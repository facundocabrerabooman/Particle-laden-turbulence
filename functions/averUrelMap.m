function [averMatrix, stdMatrix,countMatrix] = averUrelMap(allMatrices)
    % averMeanFields - Compute the average matrix from a cell array of 3D matrices.
    %
    % Syntax:
    %   averageMatrix = computeAverageMatrix(allMatrices)
    %
    % Inputs:
    %   allMatrices - Cell array containing 3D matrices (each matrix can have NaN values).
    %
    % Output:
    %   averageMatrix - 3D matrix representing the average value of each element
    %                   over all input matrices, ignoring NaN values.
    %

    % Convert the cell array of matrices into a 4D array
    stackedMatrices = cat(4, allMatrices{:});
    sz = size(stackedMatrices);
    countMatrix = zeros(sz(1),sz(2));
    % Compute the average value of each element over all matrices, ignoring NaN
    averMatrix = nanmean(stackedMatrices, [4 3]);
    stdMatrix = nanstd(stackedMatrices,0, [4 3]);
    for i = 1:sz(1)
        for j = 1:sz(2)
            countMatrix(i,j) = numel(stackedMatrices(~isnan(stackedMatrices(i,j,:,:))));
        end
    end
    countMatrix(countMatrix==0)= NaN;
    
end
