function track1 = mergeStruct(track1,track2)


Fnames=fieldnames(track2);

N1 = numel(track1);
N2 = numel(track2);

for k=1:numel(Fnames)
  for k2=1:N2
      track1(N1+k2).(Fnames{k}) = track2(k2).(Fnames{k});
  end
end