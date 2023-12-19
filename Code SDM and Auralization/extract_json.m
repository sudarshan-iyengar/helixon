function [DOA_inter, P_inter] = extract_json(ids)
    % Define the file base name
    fileBaseName = 'Interpolated_points/PC%d_interpolated_res210cm.json';

    % Preallocate cell arrays to hold the data
    DOA_inter = cell(length(ids), 1); % Cell array for DOA, each row will hold data from one file
    P_inter = cell(length(ids), 1);   % Cell array for P, each row will hold data from one file

       % Loop over the different file identifiers
    for i = 1:length(ids)
        % Create the filename
        fileName = sprintf(fileBaseName, ids(i));

        % Read and decode the JSON file
        interpolatedData = jsondecode(fileread(fileName));

        % Access position and mass data and store them in the cell arrays
        DOA_inter{i} = interpolatedData.mass;
        P_inter{i} = interpolatedData.pos;
    end
end