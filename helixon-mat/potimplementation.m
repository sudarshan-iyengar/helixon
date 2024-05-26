%% section for data handling

resolutions = [2, 3, 4];
mu = 0.01;
base_dir = 'SRIR40k/SRIR40k/';
pressure_prefix = 'ground_';
pressure_suffix = '_0_p.csv';
doa_suffix = '_0_doa.csv';
max_position = 28;

for resolution = resolutions
    for pos = 1:resolution:max_position
        if pos + resolution > max_position
            break;
        end
        
        % Load pressure and DOA data for chosen positions
        pressure_files = {sprintf('%s%s%d%s', base_dir, pressure_prefix, pos, pressure_suffix), ...
                          sprintf('%s%s%d%s', base_dir, pressure_prefix, pos + resolution, pressure_suffix)};
        doa_files = {sprintf('%s%s%d%s', base_dir, pressure_prefix, pos, doa_suffix), ...
                     sprintf('%s%s%d%s', base_dir, pressure_prefix, pos + resolution, doa_suffix)};
        
        [p_original, doa_data, p_w_pos, p_w_neg] = loadData(pressure_files, doa_files);

        % Truncate to 1000 samples for processing
        [doa_data, p_original, p_w_pos, p_w_neg] = truncateData(doa_data, p_original, p_w_pos, p_w_neg);

        % Create point clouds for positive and negative pressures
        PC1_pos = createPointCloud(doa_data{1}, p_w_pos{1});
        PC2_pos = createPointCloud(doa_data{2}, p_w_pos{2});
        PC1_neg = createPointCloud(doa_data{1}, p_w_neg{1});
        PC2_neg = createPointCloud(doa_data{2}, p_w_neg{2});
        PC1_original = createPointCloud(doa_data{1}, p_original{1});
        PC2_original = createPointCloud(doa_data{2}, p_original{2});

        % Calculate cost matrix
        C = calculateCostMatrix(PC1_pos, PC2_pos, PC1_original, PC2_original, mu);

        % Partial Optimal Transport
        distEx = 0.08 * resolution; % 8 cm per resolution as per measurements
        distTol = 0.01; % 1% tolerance chosen
        pot_pos = VSOT(PC1_pos, PC2_pos, C, distEx, distTol);
        pot_neg = VSOT(PC1_neg, PC2_neg, C, distEx, distTol);

        % Interpolate for the in-between positions
        for k = 1:resolution-1
            interp_factor = k / resolution;
            PCk_pos_better = pot_pos.interpPC(interp_factor);
            PCk_neg_better = pot_neg.interpPC(interp_factor);

            % Save interpolated results
            interp_pos = pos + k;
            saveInterpolatedResults(PCk_pos_better, PCk_neg_better, mu, interp_pos, pos, pos + resolution);
        end
        clearvars -except resolutions mu base_dir pressure_prefix pressure_suffix doa_suffix max_position resolution pos
    end
end


function [p_original, doa_data, p_w_pos, p_w_neg] = loadData(pressure_files, doa_files)
    p_original = {};
    doa_data = {};
    p_w_pos = {};
    p_w_neg = {};
    for i = 1:length(pressure_files)
        df_pressure = readmatrix(pressure_files{i});
        df_pressure = df_pressure(401:20000); % Adjust indices according to your data
        disp("PRESSURE SHAPE: ");
        disp(size(df_pressure));
        disp(df_pressure(1));
        p_original{i} = df_pressure;
        
        df_doa = readmatrix(doa_files{i});
        df_doa = df_doa(401:20000,:); % Adjust indices according to your data
        disp("DOA SHAPE: ");
        disp(size(df_doa));
        disp(df_doa(1,:));

        nan_mask = any(isnan(df_doa), 2);
        df_doa(nan_mask, :) = [];
        df_pressure(nan_mask) = [];
        
        p_original{i} = df_pressure;
        doa_data{i} = df_doa;
        
        % Splitting positive and negative pressure
        [p_pos, p_neg] = split_pos_neg(df_pressure);
        p_w_pos{i} = p_pos;
        p_w_neg{i} = p_neg;
        doa_data{i} = df_doa;
    end
end

function [doa_data, p_original, p_w_pos, p_w_neg] = truncateData(doa_data, p_original, p_w_pos, p_w_neg)
    for i = 1:2
        doa_data{i} = doa_data{i}(1:1000,:);
        p_original{i} = p_original{i}(1:1000);
        p_w_pos{i} = p_w_pos{i}(1:1000);
        p_w_neg{i} = p_w_neg{i}(1:1000);
    end
end

function PC = createPointCloud(pos, mass)
    PC = struct('pos', pos, 'mass', mass, 'n', length(mass));
end

function C = calculateCostMatrix(PC1_pos, PC2_pos, PC1_original, PC2_original, mu)
    C_euclidean = pdist2(PC1_pos.pos, PC2_pos.pos, "squaredeuclidean");
    sign_mask = (PC1_original.mass .* PC2_original.mass) < 0;
    pressure_diff = abs(PC1_original.mass - PC2_original.mass);
    penalty = mu * pressure_diff .* sign_mask;
    C = C_euclidean + penalty;
end

function saveInterpolatedResults(PCk_pos_better, PCk_neg_better, mu, interp_pos, start_pos, end_pos)
    

    filename_p_pos = sprintf('auralization/potMatMu0100_%d_%d-%d_p_pos.csv', interp_pos, start_pos, end_pos);
    filename_doa_pos = sprintf('auralization/potMatMu0100_%d_%d-%d_doa_pos.csv', interp_pos, start_pos, end_pos);
    filename_p_neg = sprintf('auralization/potMatMu0100_%d_%d-%d_p_neg.csv', interp_pos, start_pos, end_pos);
    filename_doa_neg = sprintf('auralization/potMatMu0100_%d_%d-%d_doa_neg.csv', interp_pos, start_pos, end_pos);

    writematrix(PCk_pos_better.mass, filename_p_pos);
    writematrix(PCk_pos_better.pos, filename_doa_pos);
    writematrix(PCk_neg_better.mass, filename_p_neg);
    writematrix(PCk_neg_better.pos, filename_doa_neg);
end

function visualizeIR(PCk_pos_better, PCk_neg_better)
    all_doas = [PCk_pos_better.pos; PCk_neg_better.pos];
    all_pressures = [PCk_pos_better.mass; (PCk_neg_better.mass) * -1];

    c = 343;
    timeshifts = vecnorm(all_doas, 2, 2) / c;
    time = linspace(0, 1, 100000);
    ir = zeros(size(time));

    for i = 1:length(all_pressures)
        index = round(timeshifts(i) * length(time)) + 1;
        if index <= length(ir)
            ir(index) = ir(index) + all_pressures(i);
        end
    end

    figure;
    plot(time, ir);
    title("IR");
    grid on;
end

% Function to split positive and negative pressures
function [p_pos_full, p_neg_full] = split_pos_neg(p)
    p_pos_full = zeros(size(p));
    p_neg_full = zeros(size(p));

    p_pos_indices = find(p > 0);
    p_neg_indices = find(p < 0);

    p_pos_full(p_pos_indices) = p(p_pos_indices);
    p_neg_full(p_neg_indices) = p(p_neg_indices);

    p_neg_full = abs(p_neg_full);
end
