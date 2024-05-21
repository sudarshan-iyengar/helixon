%% Define the suffixes of the impulse responses to load
baseFolder = 'Aur_Team1_noInterp30k/'; 
suffixes = { '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27'};
%% Initialize a cell array to hold the impulse responses
% Each cell will contain the IRs from all mics for a particular suffix
auralize_locs = cell(length(suffixes), 1);
fs = []; % Initialize sampling frequency

%% Load the impulse responses for each suffix and microphone
for i = 1:length(suffixes)
    auralize_set = []; % Initialize an empty array for this set of IRs
    %for j = 1:2
        % Construct the file name
        fileName = ['Auralized_SRIR_Location_Interpolated', suffixes{i}, '.wav'];
        % Construct the full path to the impulse response file
        fullPath = [baseFolder, fileName];
        % Check if the file exists
        if exist(fullPath, 'file')
            % Read the impulse response file
            [tempAu, tempFS] = audioread(fullPath);
            % Check if fs is consistent across files
            if isempty(fs)
                fs = tempFS;
            elseif fs ~= tempFS
                error(['Sampling frequency mismatch in file: ', fullPath]);
            end
            % Append to the set of IRs for this suffix
            auralize_set = [auralize_set, tempAu]; % Concatenate horizontally
        else
            warning(['File does not exist: ', fullPath]);
        end
    % end
    % Store the set of IRs in the cell array
    auralize_locs{i} = auralize_set;
end
%% For-loop to combine the files 
winlen = 33060/2;
output = auralize_locs{1};
for idx = 2:length(auralize_locs)
    output = equalize_norm(output,auralize_locs{idx},winlen);
end

savename = sprintf('int.wav');
    
    disp('Started writing the auralized wav file')
    disp([savename  ' on the disk.']);tic
    audiowrite(savename,output,48000) % <---- save the result as wav ; /10 is removed
    info = audioinfo(savename);
    disp('Wrote ... ');
    disp(info)
    disp(['... in ' num2str(toc) ' seconds'])
%%
