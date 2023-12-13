function [traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout,tracklong,Ine)

%%% Split Drop Tower Data

%load([folderin filesep 'output_post_processing.mat'],'tracklong','Ine')
%load([folderin filesep 'output_post_processing.mat'])
tracklong = tracklong(Ine);

idx_fullg = [];
idx_dec = []; % deceleration
idx_ddt = [];
% vel_fullg = [];
% vel_dec = [];
% vel_ddt = [];
% acc_fullg = [];
% acc_dec = [];
% acc_ddt = [];
traj_fullg = [];
traj_ddt = [];
traj_dec = [];

for ii = 1:numel(tracklong)
ii/numel(tracklong)
    %vel(ii) = mean([mean(tracklong(ii).Vz) mean(tracklong(ii).Vy) mean(tracklong(ii).Vx)]);
    %vel = [vel;tracklong(ii).Vz];
    %t = [t;tracklong(ii).Tf];

    idx_fullg = find(tracklong(ii).Tf<(7.29-4.1)); % 2.1 of drop + 2s after drop
    %     %full g
    if ~(isempty(idx_fullg))
        traj_fullg = [traj_fullg tracklong(ii)];
        % vel_fullg = [vel_fullg;tracklong(ii).Vz(idx_fullg)];
        % acc_fullg = [acc_fullg;tracklong(ii).Az(idx_fullg)];
    elseif ~isempty(find(tracklong(ii).Tf>(7.29-4.1) & tracklong(ii).Tf<(7.29-2)))
        traj_ddt = [traj_ddt tracklong(ii)];
        %idx_ddt = find(tracklong(ii).Tf>(7.29-4.1) & tracklong(ii).Tf<(7.29-2));
        %  vel_ddt = [vel_ddt;tracklong(ii).Vz(idx_ddt)];
        %   acc_ddt = [acc_ddt;tracklong(ii).Az(idx_ddt)];
    else
        traj_dec = [traj_dec tracklong(ii)];
        %idx_dec = find(tracklong(ii).Tf>(7.29-2));
        %vel_dec = [vel_dec; tracklong(ii).Vz(idx_dec)];
        %acc_dec = [acc_dec; tracklong(ii).Az(idx_dec)];
    end


end
clear ii idx_dec idx_ddt idx_fullg


%save([folderout filesep 'output_post_processing.mat'],'traj_dec','traj_ddt','traj_fullg','-append')
