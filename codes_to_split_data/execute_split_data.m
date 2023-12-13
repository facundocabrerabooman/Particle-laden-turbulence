%% Split Drop Tower data 
folderin = '/Volumes/landau1/Tracers/';
folderout = folderin;

% load('trajs_july3.mat')
% 
% [traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout,tracklong,Ine);
% 
% save(['dec_' fname],'traj_dec')
% save(['fullg_' fname],'traj_fullg')
% save(['ddt_' fname],'traj_ddt')

fname = 'july7a';
load('trajs_july7a.mat')

[traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout,tracklong,Ine);

save(['dec_' fname],'traj_dec')
save(['fullg_' fname],'traj_fullg')
save(['ddt_' fname],'traj_ddt')

fname = 'july7b';
load('trajs_july7b.mat')

[traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout,tracklong,Ine);

save(['dec_' fname],'traj_dec')
save(['fullg_' fname],'traj_fullg')
save(['ddt_' fname],'traj_ddt')

%%
folderin = '/Volumes/landau1/Tracers/';
folderout = folderin;

fname = 'july7c';
load('trajs_july7c.mat')

[traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout,tracklong,Ine);

save(['dec_' fname],'traj_dec')
save(['fullg_' fname],'traj_fullg')
save(['ddt_' fname],'traj_ddt')

fname = 'july9b';
load('trajs_july9b.mat')

[traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout,tracklong,Ine);

save(['dec_' fname],'traj_dec')
save(['fullg_' fname],'traj_fullg')
save(['ddt_' fname],'traj_ddt')

fname = 'july3';
load('trajs_july3.mat')

[traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout,tracklong,Ine);

save(['dec_' fname],'traj_dec')
save(['fullg_' fname],'traj_fullg')
save(['ddt_' fname],'traj_ddt')


