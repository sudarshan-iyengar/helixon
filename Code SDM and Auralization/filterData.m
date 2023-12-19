function [filtered_DOA, filtered_P] = filterData(DOA, P, threshold)
    % Filter out the samples in P below the threshold
    validIndices = P >= threshold;
    
    % Apply the valid indices to filter both P and DOA
    filtered_P = P(validIndices);
    filtered_DOA = DOA(validIndices, :); % Keep all columns (:)
end