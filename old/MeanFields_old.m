clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/exports/fullg/';

folderout = '/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/fullg/';
mkdir(folderout)
cd(folderout)

Fs=2996; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Concatenate data

if pi==pi
    trajs_conc = [];
%%% if fullg
    load([folderin filesep 'trajsf_Tracers_01_fullg.mat'],'tracklong')
    trajs_conc = [trajs_conc tracklong];

    load([folderin filesep 'trajsf_Tracers_04_fullg.mat'],'tracklong')
    trajs_conc = [trajs_conc tracklong];

%%% if ddt

    %  load([folderin filesep 'trajsf_Tracers_01_ddt.mat'],'tracklong')
    %   trajs_conc = [trajs_conc tracklong];
    % 
    %   load([folderin filesep 'trajsf_Tracers_02_ddt.mat'],'tracklong')
    %   trajs_conc = [trajs_conc tracklong];
    % 1
    % 
    % load([folderin filesep 'trajsf_Tracers_03_ddt.mat'],'tracklong')
    % trajs_conc = [trajs_conc tracklong];
    % 
    % load([folderin filesep 'trajsf_Tracers_04_ddt.mat'],'tracklong')
    % trajs_conc = [trajs_conc tracklong];
    2
2
    Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);

end
trajs_conc = trajs_conc(Ine);
clear tracklong traj_fullg

%%% change reference system
trajs_conc_tmp = trajs_conc;
counter=0;
for i=1:numel(trajs_conc)

        trajs_conc(i).x = trajs_conc_tmp(i).x;
        trajs_conc(i).y = trajs_conc_tmp(i).z;
        trajs_conc(i).z = -trajs_conc_tmp(i).y;

        trajs_conc(i).Xf = trajs_conc_tmp(i).Xf;
        trajs_conc(i).Yf = trajs_conc_tmp(i).Zf;
        trajs_conc(i).Zf = -trajs_conc_tmp(i).Yf;

        trajs_conc(i).Vx = trajs_conc_tmp(i).Vx;
        trajs_conc(i).Vy = trajs_conc_tmp(i).Vz;
        trajs_conc(i).Vz = -trajs_conc_tmp(i).Vy;

        trajs_conc(i).Ax = trajs_conc_tmp(i).Ax;
        trajs_conc(i).Ay = trajs_conc_tmp(i).Az;
        trajs_conc(i).Az = -trajs_conc_tmp(i).Ay;

        %trajs_conc(i).Z = trajs_conc(i).z;
end

%%% If doing microgravity
%disp('microgravity?'); pause

% trajs_conc_new=trajs_conc;
% for i=1:numel(trajs_conc)
%     if trajs_conc(i).t(end)<0.5
%         trajs_conc_new(i).deleteflag=1;
%     else 
%         trajs_conc_new(i).deleteflag=0;
%     end
% end
% 
% delete_index = [trajs_conc_new.deleteflag] == 1;
% trajs_conc_new(delete_index) = [];
% 
% trajs_conc=trajs_conc_new;

disp('GRAVITY goes downwards in z now; the vectors are anti-parallel')

%% Compute interpolant for MEAN flow 
dt = [4 6 8 10];
nbins = [20 20 30];
threshold = 1;
n=1;
power=1;

%[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(trajs_conc,'Xf',dt,nbins,n,power,'x','cart');
%[mYdt, mBydt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Yf',dt,nbins,n,power,'y','cart');
%[mZdt, mBzdt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Zf',dt,nbins,n,power,'z','cart');

%[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(trajs_conc,'x',dt,nbins,n,power,'x','cart');
%[mYdt, mBydt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'y',dt,nbins,n,power,'y','cart');
%[mZdt, mBzdt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'z',dt,nbins,n,power,'z','cart');
%%
[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(trajs_conc,'Vx',dt,nbins,0,1);
[mYdt, mBydt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Vy',dt,nbins,0,1);
[mZdt, mBzdt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Vz',dt,nbins,0,1);

% as dt is intergra array (like [2 4 6 8]), we have to multiply it by Fs^n
mXdt=mXdt*Fs^n;
mYdt=mYdt*Fs^n;
mZdt=mZdt*Fs^n;


% threshold = 10;
I = find(N<threshold);

mXdts=smoothn(mXdt,0.1);
mYdts=smoothn(mYdt,0.1);
mZdts=smoothn(mZdt,0.1);

mXdts(I)=NaN;
mYdts(I)=NaN;
mZdts(I)=NaN;

%%% interpolation
%[X,Y,Z]=ndgrid(bins{1},bins{2},bins{3});
[X,Y,Z]=ndgrid(linspace(-30,30,20),linspace(-30,30,20),linspace(-45,35,30));
%[X,Y,Z]=ndgrid(linspace(-100,100,50),linspace(-100,100,50),linspace(-100,100,50));

Fx=griddedInterpolant(X,Y,Z,mXdts,'cubic','none');
Fy=griddedInterpolant(X,Y,Z,mYdts,'cubic','none');
Fz=griddedInterpolant(X,Y,Z,mZdts,'cubic','none');

%save([folderout filesep 'MeanFlow_Interpolant'],'Fx','Fy','Fz')

%%% better way to plot center plane?
%figure
figure(101);clf

name = 'centerplane_betterway_xcst=0';

[y, z] = meshgrid(linspace(-30, 30, 20), linspace(-45, 35, 30));
x = 0.*ones(size(y)); % assuming the center plane is at y=0


% Evaluate interpolants on the grid
vx = Fx(x, y, z);
vy = Fy(x, y, z);
vz = Fz(x, y, z);

% Compute magnitude of velocity vector
v_mag = sqrt(vx.^2 + vy.^2 + vz.^2);

% Plot velocity vectors with colormap
quiver3(y, z, 0.*ones(size(vz)), vy, vz, zeros(size(vz)), 1,'Color', 'k'); % plot velocity vectors
hold on;
surf(y,z,v_mag, EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
cb = colorbar ; colormap jet; 
axis equal
xlabel('x(mm)')
ylabel('z(mm)')
ylabel(cb,'$|\mathbf{v}|$ (mm/s)','Interpreter','latex','Rotation',270,'FontSize',20)
view(0,90)
box; grid off

%clim([0 200]); 
%xlim([-20 20])
%ylim([-35 20])

stop
mkdir([folderout filesep 'meanfields_center'])
savefig_FC([folderout filesep 'meanfields_center' filesep name],8,6,'fig')
savefig_FC([folderout filesep 'meanfields_center' filesep name],8,6,'pdf')
stop
%ylim([-25 25])
%xlim([-15 15])

%% Compute interpolant for mean STD flow 
dt = [4 6 8 10];
nbins = [20 20 30];
threshold = 1;
n=1;
power=2;

[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(trajs_conc,'x',dt,nbins,n,power,'x','cart');
[mYdt, mBydt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'y',dt,nbins,n,power,'y','cart');
[mZdt, mBzdt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'z',dt,nbins,n,power,'z','cart');

mXdt=mXdt*Fs^n;
mYdt=mYdt*Fs^n;
mZdt=mZdt*Fs^n;

I = find(N<threshold);

mXdts=smoothn(mXdt,0.1);
mYdts=smoothn(mYdt,0.1);
mZdts=smoothn(mZdt,0.1);

mXdts(I)=NaN;
mYdts(I)=NaN;
mZdts(I)=NaN;

%%% interpolation
%[X,Y,Z]=ndgrid(bins{1},bins{2},bins{3});
[X,Y,Z]=ndgrid(linspace(-30,30,20),linspace(-30,30,20),linspace(-45,35,30));
%[X,Y,Z]=ndgrid(linspace(-100,100,50),linspace(-100,100,50),linspace(-100,100,50));

Fx=griddedInterpolant(X,Y,Z,mXdts,'cubic','none');
Fy=griddedInterpolant(X,Y,Z,mYdts,'cubic','none');
Fz=griddedInterpolant(X,Y,Z,mZdts,'cubic','none');

%save([folderout filesep 'MeanFlow_Interpolant'],'Fx','Fy','Fz')

%%% better way to plot center plane?
%figure
figure(102);clf

name = 'centerplane_xcst=0';

[y, z] = meshgrid(linspace(-30, 30, 20), linspace(-45, 35, 30));
x = 0.*ones(size(y)); % assuming the center plane is at y=0

% Evaluate interpolants on the grid
vx_std = Fx(x, y, z);
vy_std = Fy(x, y, z);
vz_std = Fz(x, y, z);

% Compute magnitude of velocity vector
v_mag_rms = sqrt(vx_std.^2 + vy_std.^2 + vz_std.^2);

% Plot velocity vectors with colormap
%quiver3(y, z, 0.*ones(size(vz)), vy, vz, zeros(size(vz)), 1,'Color', 'k'); % plot velocity vectors
hold on;
surf(y,z,v_mag_rms, EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
cb = colorbar ; colormap jet; clim([0 320]); 
axis equal
xlabel('x(mm)')
ylabel('z(mm)')
ylabel(cb,'$|\mathbf{rms(v)}|$ (mm/s)','Interpreter','latex','Rotation',270,'FontSize',20)
view(0,90)
box; grid off

xlim([-20 20])
ylim([-35 20])

mkdir([folderout filesep 'meanfields_rms_center'])
savefig_FC([folderout filesep 'meanfields_rms_center' filesep name],8,6,'fig')
savefig_FC([folderout filesep 'meanfields_rms_center' filesep name],8,6,'pdf')
stop
%ylim([-25 25])
%xlim([-15 15])


























































%% Tracer: compute the fluctuation of fluid velocity sigma_u = sqrt(<u'^2>)
dt = [4 6 8 10];
nbins = [20 20 30];
threshold = 1;

gridRange.x = [-30 30];
gridRange.y = [-30 30];
gridRange.z = [-45 35];

[gridsV,meanFieldsV,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,1,1,gridRange,1);
 [gridsVrms,meanFieldsVrms,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,1,2,gridRange,1);
 [gridsA,meanFieldsA,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,2,1,gridRange,1);
 [gridsArms,meanFieldsArms,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,2,2,gridRange,1);

sigmau1.x = sqrt(mean(vertcat(trajs_conc.sVx).^2,"omitnan"))
sigmau1.y = sqrt(mean(vertcat(trajs_conc.sVy).^2,"omitnan"))
sigmau1.z = sqrt(mean(vertcat(trajs_conc.sVz).^2,"omitnan"))
sigmau1.all = sqrt((sigmau1.x^2+sigmau1.y^2+sigmau1.z^2)/3) % unit: mm/s


%% Plot Mean Velocity -- xx->yy in meanFields
close all; clc

%plane = 12;% CENTER
plane = 2;% OFF-CENTER

mkdir([folderout filesep 'meanfields_offcenter'])

normvel = sqrt((meanFieldsV.x).^2 + (meanFieldsV.y).^2 + (meanFieldsV.z).^2);

%%% Plot center plane at constant z
name = 'vertical_plane_cst_z';
figure;

vel = reshape(normvel(:,:,plane),[25,50]);
Ys = reshape(gridsV.YY(:,:,plane),[25,50]);
Xs = reshape(gridsV.XX(:,:,plane),[25,50]);
surf(Ys,-Xs,real(vel), EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
q=quiver(gridsV.YY(:,:,plane),-gridsV.XX(:,:,plane),meanFieldsV.y(:,:,plane),-meanFieldsV.x(:,:,plane),'k','LineWidth',1);

%title(name,'Interpreter','latex')
cb = colorbar ; colormap jet; clim([0 150]); 
axis equal
%ylim([-30,30])
%xlim([-15,15])
xlabel('x(mm)')
ylabel('y(mm)')
ylabel(cb,'$|\mathbf{v}|$ (mm/s)','Interpreter','latex','Rotation',270,'FontSize',20)
view(0,90)
box; grid off

savefig_FC([folderout filesep 'meanfields_offcenter' filesep name],8,6,'fig')
savefig_FC([folderout filesep 'meanfields_offcenter' filesep name],8,6,'pdf')


%% Plot RMS Velocity
close all; clc

normvel = sqrt((meanFieldsVrms.x).^2 + (meanFieldsVrms.y).^2 + (meanFieldsVrms.z).^2);

%%% Plot center plane at constant z
name = 'rms_vertical_plane_cst_z';
figure;

vel = reshape(normvel(:,:,12),[25,50]);
Ys = reshape(gridsVrms.YY(:,:,12),[25,50]);
Xs = reshape(gridsVrms.XX(:,:,12),[25,50]);
surf(Ys,-Xs,real(vel), EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
%q=quiver(gridsVrms.YY(:,:,12),-gridsVrms.XX(:,:,12),meanFieldsVrms.y(:,:,12),-meanFieldsVrms.x(:,:,12),'k','LineWidth',1);

%title(name,'Interpreter','latex')
cb = colorbar ; colormap jet; 
clim([220 280]); 
axis equal
ylim([-30,30])
xlim([-15,15])
xlabel('x(mm)')
ylabel('y(mm)')
ylabel(cb,'$\sigma_v^2$ (mm/s)','Interpreter','latex','Rotation',270,'FontSize',20)
view(0,90)
box; grid off
stop
savefig_FC(['meanfields' filesep name],8,6,'fig')
savefig_FC(['meanfields' filesep name],8,6,'pdf')

%% Plot mean Acceleration

close all; clc

normvel = sqrt((meanFieldsA.x).^2 + (meanFieldsA.y).^2 + (meanFieldsA.z).^2);

%%% Plot center plane at constant z
name = 'mean_acc_vertical_plane_cst_z';
figure;

vel = reshape(normvel(:,:,12),[25,50]);
Ys = reshape(gridsA.YY(:,:,12),[25,50]);
Xs = reshape(gridsA.XX(:,:,12),[25,50]);
surf(Ys,-Xs,real(vel), EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
q=quiver(gridsA.YY(:,:,12),-gridsA.XX(:,:,12),meanFieldsA.y(:,:,12),-meanFieldsA.x(:,:,12),'k','LineWidth',1);

%title(name,'Interpreter','latex')
cb = colorbar ; colormap jet; 
clim([0 1e3]); 
axis equal
ylim([-30,30])
xlim([-15,15])
xlabel('x(mm)')
ylabel('y(mm)')
ylabel(cb,'$|\mathbf{a}|$ (mm/s$^2$)','Interpreter','latex','Rotation',270,'FontSize',20)
view(0,90)
box; grid off

savefig_FC(['meanfields' filesep name],8,6,'fig')
savefig_FC(['meanfields' filesep name],8,6,'pdf')

%% Plot rms Acceleration

close all; clc

normvel = sqrt((meanFieldsArms.x).^2 + (meanFieldsArms.y).^2 + (meanFieldsArms.z).^2);

%%% Plot center plane at constant z
name = 'rms_acc_vertical_plane_cst_z';
figure;

vel = reshape(normvel(:,:,12),[25,50]);
Ys = reshape(gridsArms.YY(:,:,12),[25,50]);
Xs = reshape(gridsArms.XX(:,:,12),[25,50]);
surf(Ys,-Xs,real(vel), EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
q=quiver(gridsArms.YY(:,:,12),-gridsArms.XX(:,:,12),meanFieldsArms.y(:,:,12),-meanFieldsArms.x(:,:,12),'k','LineWidth',1);

%title(name,'Interpreter','latex')
cb = colorbar ; colormap jet; 
%clim([0 1e3]); 
axis equal
ylim([-30,30])
xlim([-15,15])
xlabel('x(mm)')
ylabel('y(mm)')
ylabel(cb,'$\sigma_a^2$ (mm/s$^2$)','Interpreter','latex','Rotation',270,'FontSize',20)
view(0,90)
box; grid off

savefig_FC(['meanfields' filesep name],8,6,'fig')
savefig_FC(['meanfields' filesep name],8,6,'pdf')

