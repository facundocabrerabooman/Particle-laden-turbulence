function [neighborAll, neighborLayer] = neighbor(particle_part,tracer_tracks,Rmin,Rmax,kexp)

%%
Nframemax = 1e6;

%%

part = convertTrack(particle_part);
tracer = convertTrack(tracer_tracks);

%%
idxp = find(part.Tf<kexp*Nframemax & part.Tf>(kexp-1)*Nframemax);

for i = 1:numel(idxp)
    idx1 = idxp(i);

    %% find neighbors of current frame
    idx2 = find(tracer.Tf==part.Tf(idx1));
    [neighborAll(i,:), neighborLayer(i,:)] = neighborIdx(part,tracer,idx1,idx2,Rmin,Rmax);
    
end
