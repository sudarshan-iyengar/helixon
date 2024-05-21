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

%create an array of length winLen that linearly goes from 0 to 1
zero_to_one = linspace(0, 1, winLen);
zero_to_one = [zero_to_one; zero_to_one]';
one_to_zero = linspace(1, 0, winLen);
one_to_zero = [one_to_zero; one_to_zero]';

m1(end-winLen+1:end, :) = m1_segment.*one_to_zero + m2_segment.*zero_to_one;

xm_comp = [m1;m2(winLen+1:end, :)];
% % Overlap-add parameters
% win = hanning(winLen); % Low-pass window
% win = win/max(win);    % Normalization
% nfft = winLen;
% 
% % Process each column
% for col = 1:size(m1, 2)
%     % Apply windowing and FFT
%     M1_col = fft(bsxfun(@times, win, m1_segment(:, col)), nfft*2); 
%     M2_col = fft(bsxfun(@times, win, m2_segment(:, col)), nfft*2); 
% 
%     % Compute the difference for each column
%     comp = sqrt(abs(M2_col).^2 ./ (abs(M1_col).^2 + eps));  
%     Hf = sqrt(comp).*sqrt(conj(comp)); 
% 
% 
%     % Calculate the compensated version of the column
%     Yy_col = bsxfun(@times, Hf, M1_col); 
%     tempy_col = real(ifft(Yy_col)); 
% 
%     % Replace the equalized result in the original m1
%     m1(end-winLen+1:end, col) = m1(end-winLen+1:end, col) + tempy_col(1:winLen); 
% end
% 
% % Delete the first winLen samples in m2
% m2_updated = m2(winLen+1:end, :);
% 
% % Concatenate the updated m1 and m2
% xm_comp = [m1; m2_updated];

end




% function xm_comp = equalize_norm(m1, m2, winLen)
% % Function Description Updated
% % xm_comp = equalize_norm(m1, m2, winLen)
% % Equalizes the last winLen samples of each column in m1 to the 
% % first winLen samples of the corresponding column in m2, replaces them
% % in m1, deletes them in m2, and concatenates the results.
% 
% % PARAMETERS:
% % m1: Matrix of impulse responses [R C]
% % m2: Matrix of reference impulse responses [R C]
% % winLen: Window length in samples
% [L, N] = size(m1);
% if L < winLen
%     m2(1:winLen-L, :) = [];
%     winLen = L;
% end
% % Extract the relevant segments from m1 and m2
% m1_segment = m1(end-winLen+1:end, :); % Last winLen samples of m1
% m2_segment = m2(1:winLen, :);         % First winLen samples of m2
% 
% 
% % Overlap-add parameters
% win = hanning(winLen); % Low-pass window
% win = win/max(win);    % Normalization
% nfft = winLen;
% 
% % Process each column
% for col = 1:size(m1, 2)
%     % Apply windowing and FFT
%     M1_col = fft(bsxfun(@times, win, m1_segment(:, col)), nfft*2); 
%     M2_col = fft(bsxfun(@times, win, m2_segment(:, col)), nfft*2); 
% 
%     % Compute the difference for each column
%     comp = sqrt(abs(M2_col).^2 ./ (abs(M1_col).^2 + eps));  
%     Hf = sqrt(comp).*sqrt(conj(comp)); 
% 
% 
%     % Calculate the compensated version of the column
%     Yy_col = bsxfun(@times, Hf, M1_col); 
%     tempy_col = real(ifft(Yy_col)); 
% 
%     % Replace the equalized result in the original m1
%     m1(end-winLen+1:end, col) = m1(end-winLen+1:end, col) + tempy_col(1:winLen); 
% end
% 
% % Delete the first winLen samples in m2
% m2_updated = m2(winLen+1:end, :);
% 
% % Concatenate the updated m1 and m2
% xm_comp = [m1; m2_updated];
% 
% end


