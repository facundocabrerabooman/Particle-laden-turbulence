function newtrack = convertTrack(track)

% convert tracks to 1*1 structure array

fields = fieldnames(track);
for i = 1:length(fields)
    newtrack.(fields{i}) = vertcat(track.(fields{i}));
end