classdef VSOT < handle
    %VSOT Virtual Source Optimal Transport (aka. PCOT v2.0)
    %   Finds optimal transport couplings and interpolates between point clouds
    %   v2.0 changes the cost function to be relative to expected movement!
    % Created by: Aaron Geldert
    % Last modified: 5 Dec
    
    properties
        PC1 % point cloud 1 (struct with .pos, .mass, .n members)
        PC2 % point cloud 2
        T % transport matrix
        Tx % extended with virtual masses
        C % cost matrix
        Cx % extended with dummy cost
        distEx % expected distance (m) of transported image sources
        distTol % how far away from distEx is tolerated, as a ratio 
        dumCost % this is the cost at distRx*distTol squared
        sRatOpt % optimal ratio from optimization
        sOpt % optimal amount of mass to transport
        
        PCk % last interpolated PC result
        k % last interpolated k value
        
        % Greedy strategy of Nearest Neighbor map (baseline)
        Tg
        PCg
    end
    
    methods
        function obj = VSOT(PC1, PC2, distEx, distTol)
            %VSOT Constructor
            %   Takes in two point clouds and the virtual cost parameter
            
            % copy in data
            obj.PC1 = PC1;
            obj.PC2 = PC2;
            obj.distEx = distEx;
            obj.distTol = distTol;
            
            % Prepare cost matrix (distance-based)
            q = 2;
            mu = 0.5;
            %C = zeros(PC1.n, PC2.n);
            C = pdist2(PC1.pos, PC2.pos,"squaredeuclidean");
            sign_mask = (PC1.mass(:)*PC2.mass(:)')<0;
            pressure_diff = abs(PC1.mass(:)-PC2.mass(:)');
            penalty = mu*pressure_diff*sign_mask;
            C = C+penalty;
            disp(size(C))
            % for ii = 1:PC1.n
                %for jj = 1:PC2.n
                    %C(ii,jj) = (norm(PC1.pos(ii,:) - PC2.pos(jj,:)) - obj.distEx).^q; % normal cost
                %end
            %end
            A = max(C,[],'all') + 1; % A > largest element of C
            
            % determine the dummy cost 
            obj.dumCost = (obj.distEx*obj.distTol)^2;
            
            % add virtual point costs
            Cv1 = obj.dumCost * ones(PC1.n, 1); 
            Ctemp = horzcat(C, Cv1);
            Cv2 = obj.dumCost * ones(PC2.n, 1).'; 
            Cx = vertcat(Ctemp, horzcat(Cv2, A));
            obj.Cx = Cx;
            obj.C = C;

            % sRatio optimizer
            sVec = (1:-0.01:0).';
            sLen = numel(sVec);
            costCoarse = zeros(sLen, 1);
            
            for sInd = 1:sLen
                
                % determine max amount of mass (s)
                maxS = min(sum(abs(PC1.mass)), sum(abs(PC2.mass)));
                s = sVec(sInd) * maxS;
                
                % set virtual point masses
                m1 = [PC1.mass; sum(abs(PC2.mass)) - s];
                m2 = [PC2.mass; sum(abs(PC1.mass)) - s];
            
                % Optimization
                cvx_begin quiet
                    cvx_solver SDPT3 % SDPT3 or sedumi
                    variable Tx(PC1.n+1, PC2.n+1)
                    variable Z(PC1.n+1, PC2.n+1)
                    minimize sum(Z(:))
                    subject to
                        Tx >= 0; % matrix must be nonnegative
                        sum(sum(Tx(1:PC1.n, 1:PC2.n))) == s; % total mass transported, w/o virtual pts
                        Tx * ones(PC2.n+1,1) == m1; % marginal sums kept, including virtual points
                        Tx' *ones(PC1.n+1,1) == m2;
                        Z == Cx .* Tx;
                cvx_end

                Tx(end, end) = 0; % this must be true
                Tx(Tx < 1e-5) = 0; % round down to 0
                if isnan(Tx)
                    warning('Transport matrix contains NaN elements.');
                end
                
                costCoarse(sInd) = sum(sum(Cx.*Tx));
                if (sInd>1) && (costCoarse(sInd) > costCoarse(sInd-1))
                    % stopping condition, don't update optimal values
                    disp(['Found optimal sRatio = ' num2str(obj.sRatOpt)]);
                    break;
                else
                    % improved, save data
                    obj.Tx = Tx;
                    obj.T = Tx(1:PC1.n, 1:PC2.n);
                    obj.sRatOpt = sVec(sInd);
                    obj.sOpt = obj.sRatOpt * maxS;
                end
                
            end % end for
            
            % GREEDY transport map assignment algorithm:
            % How does it handle unbalanced mass?
            % Normalization for transport map; rescale for interp
            
            unmapped1 = PC1.mass ./ sum(PC1.mass);
            unmapped2 = PC2.mass ./ sum(PC2.mass);
            [~, sortInds] = sort(obj.C(:));
            [i, j] = ind2sub(size(C),sortInds);
            obj.Tg = zeros(size(C));
            % go through sorted entries of C
            for ii = 1:length(sortInds)
                % assign as much unmapped mass as possible
                mass = min(unmapped1(i(ii)), unmapped2(j(ii)));
                obj.Tg(i(ii),j(ii)) = mass;
                unmapped1(i(ii)) = unmapped1(i(ii)) - mass;
                unmapped2(j(ii)) = unmapped2(j(ii)) - mass;
                
                % stopping criterion
                if min(sum(unmapped1),sum(unmapped2) == 0)
                    break;
                end
            end
            
        end % end function
        
        
        function PCk = interpPC(obj, k)
            %INTERPC Interpolates the partial OT problem at k
            
            assert(k>=0 && k<=1, 'k must be between 0 and 1.');
            
            %% There are 3 types of mass to account for:
            % 1. Transported Mass (TM)
            % - The interior of the transport matrix defines the masses of all
            % couplings of real points to other real points. Their corresponding
            % positions are simply linearly interpolated
            %
            % 2. Growing Mass (GM)
            % - When a real point is coupled with a virtual point and one or 
            % more other real points, the mass of the coupling is linearly 
            % combined with the other transported masses coupled with the real point
            %
            % 3. Stationary Mass (SM)
            % - When a real point is coupled only with a virtual point, the 
            % mass of the coupling spontaneously appears at position the real 
            % point posisiton, mass proportional to the interpolation parameter

            % Recall that mass is defined using ENERGY!
            
            % 1. Transported Mass (TM)
            Ttm = obj.T; % ignore virtual mass columns

            % 2. Moving Mass (MM)
            % From X->Y 
            xRat = Ttm./vecnorm(Ttm, 1, 2); % ratio of MM distribution, rowsum=1
            Tmm1 = xRat .* obj.Tx(1:end-1, end); % calculate MM values
            Tmm1(isnan(Tmm1)) = 0; % remove NaN elements

            % From Y->X
            yRat = Ttm./vecnorm(Ttm, 1, 1); % colsum = 1
            Tmm2 = obj.Tx(end, 1:end-1) .* yRat;
            Tmm2(isnan(Tmm2)) = 0;

            % Calculate transported & growing masses at interp value k
            Ttmk = Ttm + (1-k)*Tmm1 + k*Tmm2; % fixed this!
            [I, J, TM] = find(Ttmk);
            numTm = numel(TM);

            % Calculate transported mass positions at interp value k
            posTm = zeros(numTm,3);
            for ii = 1:numTm
                posTm(ii,:) = (1-k).*obj.PC1.pos(I(ii),:) + k.*obj.PC2.pos(J(ii),:);
            end

            % 3. Stationary Mass (SM)
            % Which virtual masses are coupled only to Y?
            whereY = isnan(sum(yRat, 1));
            Sm1 = k* obj.Tx(end,whereY).'; % ordered vector
            posSm1 = obj.PC2.pos(whereY,:);

            whereX = isnan(sum(xRat, 2)); % virtual masses coupled only to X
            Sm2 = (1-k)* obj.Tx(whereX,end);
            posSm2 = obj.PC1.pos(whereX,:);

            % Concatenate masses and form point cloud struct
            PCk.pos = [posTm; posSm1; posSm2];
            PCk.mass = [TM; Sm1; Sm2];
            PCk.n = length(PCk.mass);

            % object remembers the last interpolation (for plotting)
            obj.k = k;
            obj.PCk = PCk;
        end
        
        function PCg = greedyInterp(obj, k)
            
            [ig, jg, tg] = find(obj.Tg);
            numTg = numel(tg);

            % Calculate transported mass positions at interp value k
            posTg = zeros(numTg,3);
            for ii = 1:numTg
                posTg(ii,:) = (1-k).*obj.PC1.pos(ig(ii),:) + k.*obj.PC2.pos(jg(ii),:);
            end
            
            scaleFactor = (1-k)*sum(obj.PC1.mass) + k*sum(obj.PC2.mass); 
            
            PCg.pos = posTg;
            PCg.mass = tg * scaleFactor;
            PCg.n = numel(tg);
            
            obj.k = k;
            obj.PCg = PCg;
        end
    end
end

