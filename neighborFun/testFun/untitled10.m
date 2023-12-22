

%% prepare idx for searching
lenT = arrayfun(@(x) numel(x.T),tracklong);
% Compute the cumulative sum from the 1st to each element
cumulative_sum = cumsum(lenT);

Ntrack = arrayfun(@(x) unique(x.Ntrack),tracklong);

%% transform the tracks structure
tracer_tracks = tracklong;
fields2 = fieldnames(tracer_tracks);
for i = 1:length(fields2)
    tracer.(fields2{i}) = vertcat(tracer_tracks.(fields2{i}));
end

%% substrat tracks
for i =1:numel(neighbor_global)
    for j = 1:numel(neighbor_global(i).idx)
        idx0 = neighbor_global(i).idx(j);
        ntrack = tracer(idx0).Ntrack;
        
        idx1 = find(ntrack==ntrack);

        if idx1 == 1
            track0 = tracklong(1);
            idx2 = idx0;
        else
            track0 = tracklong(idx1-1);
            idx2 = idx0-cumulative_sum(idx1-1);
        end
        fieldname = fieldnames(track0);
        for k = 1:numel(fieldname)
            tracksubs.(fieldname(k)) = track0.(fieldname(k))(idx2);
        end

        





