function [trackFront, trackBack]  = neighborAllTS(particle_part,tracer_tracks,neighborAll,startTime,endTime,kexp)

%% terminal state range
Nframemax = 1e6;

tstart = mod(startTime(kexp),Nframemax);
tend = mod(endTime(kexp),Nframemax);

%% convert the tracks structure

part = convertTrack(particle_part);
tracer = convertTrack(tracer_tracks);

%%
fieldToCopy = fieldnames(tracer_tracks);

%% prepare idx for searching
idxfront = sort(vertcat(neighborAll(tstart:tend).idxfront));
idxback = sort(vertcat(neighborAll(tstart:tend).idxback));

%% 
if ~isempty(idxfront)
    trackFront = retainNeighborTracks(part,tracer,fieldToCopy,idxfront);
else
    trackFront = [];
end
if ~isempty(idxback)
    trackBack  = retainNeighborTracks(part,tracer,fieldToCopy,idxback);
else
    trackBack = [];
end
