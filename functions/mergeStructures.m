function catStruc = mergeStructures(track_struct)
    
    % Get field names
    fieldName = fieldnames(track_struct);

    n = max(size(track_struct));

    % Initialize the structure with the first element
    catStruc = track_struct(1);
    
    % Determine the length of the first structure
    
    
    % Process the second structure
    for j = 2:n
        for k = 1:numel(fieldName) 
            len1 = numel(catStruc.(fieldName{k}));
            len2 = numel(track_struct(j).(fieldName{k}));
            % Concatenate the field along the rows for each element
            catStruc.(fieldName{k})(len1+1:len1+len2, 1) = track_struct(j).(fieldName{k});
        end
    end
end
