function [maxTimeLength,startTime,endTime,Thres] = findTerminalState(part,RelError)

fields = fieldnames(part);
for i = 1:length(fields)
    part2.(fields{i}) = vertcat(part.(fields{i}));
end

Nframemax = 1e6;
% NTrackmax = 1e9;
for i = 1:numel(part2)
    nexp(i,:) = fix(part2(i).Tf/Nframemax)+1;
%     t(i,:) = mod(part2(i).Tf,Nframemax)/fps;
end
uni_exp = unique(nexp);
Nexp = numel(uni_exp);


maxTimeLength = zeros(Nexp,1);
startTime = zeros(Nexp,1);
endTime  = zeros(Nexp,1);
Thres  = zeros(Nexp,1);

for j = 1:Nexp
    idx = find(nexp==j);
    [maxTimeLength(j), startTime(j), endTime(j),Thres(j)] = findTimeRange(part2.Ay(idx), RelError);
    startTime(j) = startTime(j) +  Nframemax *(j-1);
    endTime(j)   = endTime(j)   +  Nframemax *(j-1);
end

