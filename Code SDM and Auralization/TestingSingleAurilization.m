%% SDM toolbox for analysis of GRASVI50 impulse responses. 
clear; clc; close all;

%% set current folder to the one containing this script for relative paths 
% to work properly
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear tmp;

% add base folder to path for dependencies to be available
currentFolder = '/Users/wille/IIW/Master/R&D/SDM implementation';

%% Define the base folder and microphone prefixes
baseFolder = 'measurements/';
micPrefixes = {'micA_', 'micB_', 'micC_', 'micD_', 'micE_', 'micF_'};

%% Define the suffixes of the impulse responses to load
suffixes = {'011'};

%% Initialize a cell array to hold the impulse responses
% Each cell will contain the IRs from all mics for a particular suffix
ir_locs = cell(length(suffixes), 1);
fs = []; % Initialize sampling frequency

%% Load the impulse responses for each suffix and microphone
for i = 1:length(suffixes)
    ir_set = []; % Initialize an empty array for this set of IRs
    for j = 1:length(micPrefixes)
        % Construct the file name
        fileName = [micPrefixes{j}, suffixes{i}, '.wav'];
        % Construct the full path to the impulse response file
        fullPath = [baseFolder, fileName];
        % Check if the file exists
        if exist(fullPath, 'file')
            % Read the impulse response file
            [tempIR, tempFS] = audioread(fullPath);
            % Check if fs is consistent across files
            if isempty(fs)
                fs = tempFS;
            elseif fs ~= tempFS
                error(['Sampling frequency mismatch in file: ', fullPath]);
            end
            % Append to the set of IRs for this suffix
            ir_set = [ir_set, tempIR]; % Concatenate horizontally
        else
            warning(['File does not exist: ', fullPath]);
        end
    end
    % Store the set of IRs in the cell array
    ir_locs{i} = ir_set;
end

%% Pre-processing: Load the Chirp Signal & Deconvolution
% Load the chirp signal from a WAV file
[inv_chirp_signal, chirp_fs] = audioread('inv_chirp_signal.wav');

% Check if the sampling rate of the chirp matches that of the impulse responses
if chirp_fs ~= fs
    error('Sampling frequency mismatch between the chirp signal and impulse responses.');
end

% Define the length of the chirp signal
lenChirp = length(inv_chirp_signal);
length_ir = 48000;

% Initialize a cell array to store the computed impulse responses
impulse_responses = cell(length(ir_locs), 1);

% Perform deconvolution for each location
for idx = 1:length(ir_locs)
    % Truncate the recorded signals to the length of the chirp signal
    truncated_recordings = ir_locs{idx}(1:lenChirp, :);

    % FFT of the chirp signal
    inv_chirp_fft = fft(inv_chirp_signal, lenChirp);

    % Initialize an array to store the deconvolved responses for each microphone
    deconvolved_responses = zeros(length_ir, size(truncated_recordings, 2));

    % Deconvolve for each microphone
    for mic = 1:size(truncated_recordings, 2)
        % FFT of the truncated recorded signal
        recorded_fft = fft(truncated_recordings(:, mic), lenChirp);

        % Perform deconvolution in the frequency domain
        ir_fft = recorded_fft .* inv_chirp_fft;

        % Inverse FFT to get the impulse response
        ir_time = ifft(ir_fft, lenChirp, 'symmetric');
        deconvolved_responses(:, mic) = ir_time(1:length_ir,:);
    end

    % Store the deconvolved responses (impulse responses) for this location
    irs{idx} = deconvolved_responses;
end

%% Post-processing: calculate window length
% Given values
d_max = 0.05; % You need to fill in the maximum distance between mics in meters (50mm spacer)
c = 343; % Speed of sound in m/s at room temperature

% Calculate the minimum window length in seconds
min_winLen = 2 * d_max / c;

% Convert the minimum window length to samples
% This assumes 'fs' is already set to the sampling frequency in Hz
min_winLen_samples = ceil(min_winLen * fs);

% Set the window length to be just over the minimum
min_winLen = min_winLen_samples + 1; % or any small number of samples over the minimum
recommend_winLen = 34;

%% Post-processing: calculate parallel frames
numRows = size(ir_set, 1); % Get the length of the impulse response
n = nextpow2(numRows);   % Find the next power of two
nextPowerOfTwo = (2^n); % Calculate the next power of two

%% Post-processing: creation of the SDM struct
a = createSDMStruct('DefaultArray', 'GRASVI50', 'c', c, 'fs', fs, 'winLen', recommend_winLen, 'parFrames', nextPowerOfTwo);

%% Post-processing: calculate the SDM coefficients
% Initialize a cell array to store DOA results for each location
    
    % Calculate the DOA using the SDMPar function
DOA = SDMPar(irs{1}, a);

%% Post-processing: calcualte pressure

    % Averaging across all microphones (columns)
    P_avg = mean(irs{1}, 2);
    P = irs{1};

    %% Post-processing: struct for visualization with a set of parameters
% Load default setup for a small room and change some of the variables
v = createVisualizationStruct('DefaultRoom','Small',...
    'name','Acoustic Lab','fs',fs);
% For visualization purposes, set the text interpreter to latex
set(0,'DefaultTextInterpreter','latex')


%% Pressure
% Convert each cell of P_avg into a cell array containing a single element
P_avg_modified= {P_avg};  % Convert each double array into a cell array

% Now try visualizing with the modified structure

fig = figure;  % Create a new figure for each location
parameterVisualization(P_avg_modified, v);  % Visualize using the modified data structure
title(sprintf('Location %d', idx));  % Set a title for each plot
% Save the figure
close(fig);


%% Visualizing: Draw the spatio temporal visualization for each section plane
% Iterate over each location

% Select the pressure data and DOA for the current location
current_P = P;
current_DOA = DOA;  % Similarly, extract the array from the cell

% Create spatio-temporal visualizations for the current location
% Adjust v.plane as needed for each visualization

% Lateral plane visualization
v.plane = 'lateral';
fig = figure;  % Create a new figure for the visualization
spatioTemporalVisualization({current_P}, {current_DOA}, v);
title(sprintf('Location %d - Lateral Plane'));
close(fig);

% Transverse plane visualization
v.plane = 'transverse';
fig = figure;  % Create a new figure for the visualization
spatioTemporalVisualization({current_P}, {current_DOA}, v);
title(sprintf('Location %d - Transverse Plane'));
close(fig);

% Median plane visualization
v.plane = 'median';
fig = figure;  % Create a new figure for the visualization
spatioTemporalVisualization({current_P}, {current_DOA}, v);
title(sprintf('Location %d - Median Plane'));
close(fig);


% % <----- EOF SRIR

%% Download a KU-100 HRIR from the Cologne audio team server 
% HRTF is the fourier transform of HRIR
% (skipped automatically if the HRIR dataset already exists)
HRIR_URL = 'https://zenodo.org/record/3928297/files/HRIR_FULL2DEG.sofa?download=1';
HRIR_Folder   = currentFolder;
HRIR_Filename = 'KU100_HRIR_FULL2DEG_Koeln.sofa';
HRIR_Subject  = 'KU100';    % Name of the HRIR subject (only used for naming purposes while saving).
HRIR_Type     = 'SOFA';     % File format of the HRIR. Only SOFA is supported for now.

[~,~] = mkdir(HRIR_Folder); % ignore warning if directory already exists
HRIR_Path = fullfile(HRIR_Folder, HRIR_Filename);
clear HRIR_Folder HRIR_Filename;

fprintf('Downloading HRIR dataset ... ');
if isfile(HRIR_Path)
    fprintf('skipped.\n\n');
else
    websave(HRIR_Path, HRIR_URL);
    fprintf('done.\n\n');
end; clear HRIR_URL;

%% Load HRIR and Process into HRTF from the KU100 Dummy Head
% Ensure you have SOFA_Toolbox_for_Matlab_and_Octave_2 installed
% SOFAtoolbox: https://github.com/sofacoustics/SOFAtoolbox
% Add the SOFA Toolbox to MATLAB Path:
    % Open MATLAB.
    % In the MATLAB environment, go to the "Home" tab.
    % Click on "Set Path" (this opens a dialog where you can add or remove paths).
    % Click on "Add with Subfolders...". in /Users/jeffeehsiung/Library/Application Support/MathWorks/MATLAB Add-Ons/Toolboxes/
    % Navigate to the directory where you unzipped the SOFA Toolbox.
    % Select the root folder of the SOFA Toolbox and click "OK".
    % Click "Save" in the Set Path dialog to save the changes.
    % Click "Close" to close the dialog.
% In the MATLAB Command Window, type SOFAstart and press Enter.
SOFAstart;
hrtfData = SOFAload('KU100_HRIR_FULL2DEG_Koeln.sofa');
% Extract HRTF and source position information
hrtfs = hrtfData.Data.IR;  % HRTF signals
sourcePositions = hrtfData.SourcePosition;  % Corresponding source positions
% Display some information about the impulse response
SOFAinfo(hrtfData);
% Have a look at the size of the data
disp(['size [MxRxN]: ' num2str(size(hrtfs))])


%% Load the music signal
% Read stereo signal
[S, fs] = audioread('120 BPM - ROCK.wav');

% Choose 10 seconds and resample
%Sr = resample(S(1:44.e3*10,:),480,441); % Adjust the duration and sampling rate as needed


%% Create a struct for synthesis with a set of parameters
% Create a struct for synthesis with your SRIRs
s = createSynthesisStruct('Binaural',true,...
        'DefaultArray', 'KU100', ...
        'hrtfData',hrtfData,... 
        'snfft',length(P_avg),...
        'fs',fs, ...
        'c',c, ...
        'Radius', 3.5);


%% Iterate over each location for SRIR

% Select the SRIR data for the current location
current_SRIR = P_avg;
current_DOA = DOA;

% Synthesize the spatial impulse response with KU100 HRIR
Hbin = cell(1,2);
for channel = 1:2
    [~, Hbin{channel}] = synthesizeSDMCoeffs(current_SRIR, current_DOA, s);
end
disp(['Finished synthesizing SIR']);
% Convolution with an audio signal (can be a test signal or music)
Sr = resample(S, fs, hrtfData.Data.SamplingRate);  % Make sure the audio signal is resampled to match the HRIR sampling rate
Y = zeros(size(Sr, 1), 2);
for channel = 1:2  % Measured left and right channel
    for ear = 1:2  % Left and right ear
        Y(:, ear) = Y(:, ear) + fftfilt(Hbin{channel}(:, ear), Sr(:, channel));
    end
end
disp(['Finished convolution']);
% Save the auralized signal for this location
% Save the file to the default folder with a custom filename.
% Save the result as wav, as wav can handle upto 256 channels.
disp('Started Auralization');tic

savename = sprintf('Auralized_SRIR_Location1.wav');
if max(abs(Y(:))) > 1
    Y = Y/max(abs(Y(:)))*.9;
    disp('Sound normalized, since otherwise would have clipped')
end
disp(['Ended Auralization in ' num2str(toc) ' seconds.'])
disp('Started writing the auralization wav file')
disp([savename  ' on the disk.']);tic
audiowrite(savename,Y/10,s.fs) % <---- save the result as wav
info = audioinfo(savename);
disp('Wrote ... ');
disp(info)
disp(['... in ' num2str(toc) ' seconds'])


%% test:

% Select the SRIR data for the current location
current_SRIR = P_avg;
current_DOA = DOA;

% Synthesize the spatial impulse response with KU100 HRIR
Hbin = cell(1,2);
for channel = 1:2
   % [~, Hbin{channel}] = synthesizeSDMCoeffs(current_SRIR{channel}, current_DOA, s);
    [~, Hbin{channel}] = synthesizeSDMCoeffs(current_SRIR, current_DOA, s);

end
disp(['Finished synthesizing SIR']);
% Convolution with an audio signal (can be a test signal or music)
Sr = resample(S, fs, hrtfData.Data.SamplingRate);  % Make sure the audio signal is resampled to match the HRIR sampling rate
Y = zeros(size(Sr, 1), 2);
for channel = 1:2  % Measured left and right channel
    for ear = 1:2  % Left and right ear
        Y(:, ear) = Y(:, ear) + fftfilt(Hbin{channel}(:, ear), Sr(:, channel));
    end
end
disp(['Finished convolution']);
% Save the auralized signal for this location
% Save the file to the default folder with a custom filename.
% Save the result as wav, as wav can handle upto 256 channels.
disp('Started Auralization');tic
%%
savename = sprintf('Auralized_SRIR_Location%d.wav', 123);
if max(abs(Y(:))) > 1
    Y = Y/max(abs(Y(:)))*.9;
    disp('Sound normalized, since otherwise would have clipped')
end
disp(['Ended Auralization in ' num2str(toc) ' seconds.'])
disp('Started writing the auralization wav file')
disp([savename  ' on the disk.']);tic
audiowrite(savename,Y/10,s.fs/4) % <---- save the result as wav
info = audioinfo(savename);
disp('Wrote ... ');
disp(info)
disp(['... in ' num2str(toc) ' seconds'])
