function mergedStruc = mergeStructures(struct1, struct2)
    % Merge two structures by concatenating them along fields
    
    % Get field names
    fieldNames1 = fieldnames(struct1);
    fieldNames2 = fieldnames(struct2);

    % Ensure both structures have the same fields
    assert(isequal(fieldNames1, fieldNames2), 'Structures must have the same fields.');

    % Initialize the grouped structure with the first structure
    mergedStruc = struct1;
    
    % Determine the length of the first structure
    len1 = numel(struct1);
    
    % Process the second structure
    for j = 1:numel(struct2)
        % Loop through each field in the second structure
        for i = 1:numel(fieldNames2)
            fieldName = fieldNames2{i};
            
            % Concatenate the field along the rows for each element
            mergedStruc(len1+j, 1).(fieldName) = struct2(j).(fieldName);
        end
    end
end
