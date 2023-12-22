
%% Calculate MSD
f_fps=[0.354 3125 ; 0.500 3125 ; 0.707 3125 ; 1.000 3125 ; 1.414 6250 ; 2.000 6250 ; 2.378 6250];

for k_ffps = 1 : 1%size(f_fps,1)
    prefix = ['tracked_f' num2str(0.354,'%1.3f') '_' num2str(3125) 'fps_'];  
    dd = dir([prefix '*.bin']);
    
    for kfile = 1 : numel(dd)
        part(kfile) = readTracksLEM(dd(kfile).name,kfile);
    end

tracks =  part2trackLEM(part);

% msd.frot = f_fps(k_ffps,1);
% msd.fech = f_fps(k_ffps,2);
% msd.msdx = structFunc_struct(tracks,'x',2);
% msd.msdy = structFunc_struct(tracks,'y',2);
% msd.msdz = structFunc_struct(tracks,'z',2);

end

%%
%msd.frot = f_fps(k_ffps,1);
%msd.fech = f_fps(k_ffps,2);
msd2.msdx = structFunc_struct(tracksStitched3,'x',2);
msd2.msdy = structFunc_struct(tracksStitched3,'y',2);
msd2.msdz = structFunc_struct(tracksStitched3,'z',2);