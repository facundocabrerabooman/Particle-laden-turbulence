function trackout = mergeStruct(track1,track2)


Fnames=fieldnames(part);

N1 = numel(track1);
N2 = numel(track2);

for k=1:numel(Fnames)
  trackout(1:N1).(Fnames(k)) = track1.(Fnames(k));
  trackout(N1+1:N1+N2).(Fnames(k)) = track2.(Fnames(k));
end