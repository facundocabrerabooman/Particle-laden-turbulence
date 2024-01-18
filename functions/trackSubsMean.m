function track_subsMean = trackSubsMean(track)

fieldin  = {'Tf','Xf','Yf','Zf','sVx','sVy','sVz','sAx','sAy','sAz'};
fieldout = {'Tf','Xf','Yf','Zf','Vx','Vy','Vz','Ax','Ay','Az'};

for k = 1:numel(track)
    for i =1:numel(fieldin)
        track_subsMean(k).(fieldout{i}) = track(k).(fieldin{i});
    end
end

for k = 1:numel(track)
    Inan = isnan(track_subsMean(k).Vx) | isnan(track_subsMean(k).Vy) | isnan(track_subsMean(k).Az) | isnan(track_subsMean(k).Ax) | isnan(track_subsMean(k).Ay) | isnan(track_subsMean(k).Az);
    for i =1:numel(fieldout)
        track_subsMean(k).(fieldout{i})(Inan) = [];
    end
end

N = arrayfun(@(X)(numel(X.Tf)),track_subsMean);
track_subsMean(N==0) = [];
    