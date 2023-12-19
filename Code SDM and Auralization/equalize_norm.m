function xm_comp = equalize_norm(m1, m2, winLen)
% Function Description Updated
% xm_comp = equalize_norm(m1, m2, winLen)
% Equalizes the last winLen samples of each column in m1 to the 
% first winLen samples of the corresponding column in m2, replaces them
% in m1, deletes them in m2, and concatenates the results.

% PARAMETERS:
% m1: Matrix of impulse responses [R C]
% m2: Matrix of reference impulse responses [R C]
% winLen: Window length in samples
[L, N] = size(m1);
if L < winLen
    m2(1:winLen-L, :) = [];
    winLen = L;
end
% Extract the relevant segments from m1 and m2
m1_segment = m1(end-winLen+1:end, :); % Last winLen samples of m1
m2_segment = m2(1:winLen, :);         % First winLen samples of m2

% Overlap-add parameters
win = hanning(winLen); % Low-pass window
win = win/max(win);    % Normalization
nfft = winLen;

% Process each column
for col = 1:size(m1, 2)
    % Apply windowing and FFT
    M1_col = fft(bsxfun(@times, win, m1_segment(:, col)), nfft*2); 
    M2_col = fft(bsxfun(@times, win, m2_segment(:, col)), nfft*2); 

    % Compute the difference for each column
    comp = sqrt(abs(M2_col).^2 ./ (abs(M1_col).^2 + eps));  
    Hf = sqrt(comp).*sqrt(conj(comp)); 
    

    % Calculate the compensated version of the column
    Yy_col = bsxfun(@times, Hf, M1_col); 
    tempy_col = real(ifft(Yy_col)); 

    % Replace the equalized result in the original m1
    m1(end-winLen+1:end, col) = m1(end-winLen+1:end, col) + tempy_col(1:winLen); 
end

% Delete the first winLen samples in m2
m2_updated = m2(winLen+1:end, :);

% Concatenate the updated m1 and m2
xm_comp = [m1; m2_updated];

end


% function xm_comp = equalize_norm(m1,m2,winLen)
% % Function Description Updated
% % xm_comp = equalize_norm(m1, m2, winLen)
% % Matches each column in m1 to the corresponding column in m2
% % using overlap-add in windows of size winLen.
% %
% % PARAMETERS:
% % m1: Matrix of impulse responses [R C]
% % m2: Matrix of reference impulse responses [R C]
% % winLen: Window length in samples
% % 
% 
% % Add zeros to fully construct the IR
% N = size(m1,1); % N will be the number of rows
% % concatenate the matrix with dimension: [winLen x numbers columns of m1] to x_sdm (dimension: [N x numbers columns of m1]) on the top and the bottom
% m1 = [zeros(winLen,size(m1,2));m1;zeros(winLen,size(m1,2))]; 
% % concatenate the matrix with dimension: [winLen x numbers columns of x_ref] to x_sdm (dimension: [N x numbers columns of x_ref]) on the top and the bottom
% m2 = [zeros(winLen,size(m2,2));m2;zeros(winLen,size(m2,2))];
% 
% % Overlap-add parameters
% win = hanning(winLen); % low-pass
% win = win/max(win); % normalization
% nfft = winLen;
% startx = 1;
% endx = winLen;
% nforward = winLen/2;
% 
% % Number of frames, i.e., windows
% NN = round((size(m1,1)-winLen)/nforward);
% % y_dim is [(number of rows m1 + 2 * winlen) x number of columns m1]
% y = zeros(size(m1,1)+2*nfft, size(m1,2));
% for n = 1:NN
% 
%     % Disply the current frame number
%     if(mod(n,1000)==0)
%         disp(['equalizeNLS: processing frame : ' num2str(n)])
%     end
% 
%     for col = 1:size(m1, 2)
%         % Process each column separately
%         M2_col = fft(bsxfun(@times, win, m2(startx:endx, col)), nfft*2);
%         M1_col = fft(bsxfun(@times, win, m1(startx:endx, col)), nfft*2);
% 
%         % Compute the difference for each column
%         comp = sqrt(abs(M2_col).^2 ./ (abs(M1_col).^2 + eps));
%         Hf = sqrt(comp).*sqrt(conj(comp));
% 
%         % Calculate the compensated version of the column
%         Yy_col = bsxfun(@times, Hf, M1_col);
%         tempy_col = real(ifft(Yy_col));
% 
%         % Overlap-add for each column
%         y(startx:startx + nfft*2-1, col) = y(startx:startx + nfft*2-1, col) + tempy_col;
%     end
% 
%     startx = startx + nforward;
%     endx = endx + nforward;
% 
% end
% % Truncate to the original size
% xm_comp = y(winLen+1:winLen+N,:);
% 
% end
