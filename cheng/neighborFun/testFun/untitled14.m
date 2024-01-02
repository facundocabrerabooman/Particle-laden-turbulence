idxfront = sort(vertcat(neighbor_global.idxfront));
idxback = sort(vertcat(neighbor_global.idxback));

fieldname = fieldnames(tracer)
tt = struct();
for i = 1:numel(idxfront)
    for j = 1:numel(fieldname)
        tt.(fieldname{j})(i,1) = tracer.(fieldname{j})(idxfront(i));
    end
end

ntrack = unique(tt.Ntrack);
fieldname = fieldnames(tracklong)
for i = 1:numel(ntrack)
    idx = find(tt.Ntrack == ntrack(i));
    for j = 1:numel(fieldname)
        trackstr(i).(fieldname{j}) = tt.(fieldname{j})(idx);
    end
end