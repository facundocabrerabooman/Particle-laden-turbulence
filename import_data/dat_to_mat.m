%% STB .DAT to MATLAB

function d = dat_to_mat(pathin, fname)

[zone,VARlist] = tec2mat([pathin filesep fname '.dat']);
%save(data_struct)

%load([pathin filesep fname '.mat'])

%%% Transform data to use tracj3d_fc.m (i.e. get trajs)

idx = [];x=[];y=[];z=[];
Vx=[];Vy=[];Vz=[];
Ax=[];Ay=[];Az=[];
for ii=1:size(zone,2)

    ii/size(zone,2)*100
   
    idx = [idx ; ii.*ones(numel(zone(ii).data(:,1)),1)]; % index
    x = [x ; zone(ii).data(:,1)]; % x
    y = [y ; zone(ii).data(:,2)]; % y
    z = [z ; zone(ii).data(:,3)]; % z

    Vx = [Vx ; zone(ii).data(:,5)]; % x
    Vy = [Vy ; zone(ii).data(:,6)]; % y
    Vz = [Vz ; zone(ii).data(:,7)]; % z

    Ax = [Ax ; zone(ii).data(:,10)]; % x
    Ay = [Ay ; zone(ii).data(:,11)]; % y
    Az = [Az ; zone(ii).data(:,12)]; % z

end

d = [idx x y z Vx Vy Vz Ax Ay Az];

end
