function track = part2track(part)

% track = part2track(part)
%
% converts a part structure (where particles positions are stored in time
% frame order) into a track structure (with particles stored in
% trajectories order)

Fnames=fieldnames(part);
%I=find(strcmp(Fnames,'T'));
%Fnames(I)=[];

for k=1:numel(Fnames)
    if size(part(1).(Fnames{k}),2)==1
        F{k}=vertcat(part.(Fnames{k}))';
    else
        F{k}=[part.(Fnames{k})];
    end
end

if ~isempty(find(strcmp(Fnames,'Ntrack')))
    if size(part(1).Ntrack,2)==1
        Ntrack=real(vertcat(part.Ntrack)');
    else
        Ntrack=real([part.Ntrack]);
    end
elseif ~isempty(find(strcmp(Fnames,'dNtrack')))
    if size(part(1).dNtrack,2)==1
        Ntrack=(vertcat(part.dNtrack)');
    else
        Ntrack=([part.dNtrack]);
    end
else
    return
end
[Ns,Is]=sort(Ntrack);
%figure;plot(Is)
%size(Ns)

for k=1:numel(Fnames)
    size(F{k})
    F{k}=F{k}(Is);
end

Tracks = bwconncomp(1-abs(diff(Ns))>=1);
I0=cellfun(@(X)(X(1)),Tracks.PixelIdxList);
Ltrack=cellfun(@numel,Tracks.PixelIdxList);
[NtrackU]=unique(Ntrack(Is));

for ktraj=1:Tracks.NumObjects
    for k=1:numel(Fnames)
        track(ktraj).(Fnames{k})=F{k}(I0(ktraj):I0(ktraj)+Ltrack(ktraj))';
    end
end