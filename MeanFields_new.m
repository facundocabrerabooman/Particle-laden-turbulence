clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Users/fcb/Downloads/tracers/ddt/';

folderout = '/Users/fcb/Downloads/analysis_tracers/ddt/';
mkdir(folderout)
cd(folderout)

Fs=2990; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Concatenate data

if pi==pi
    trajs_conc = [];

%    load([folderin filesep 'trajsf_TrCer_1000_01_ddt_tracer.mat'],'tracklong')
%    trajs_conc = [trajs_conc tracklong];
% 1
%     load([folderin filesep 'trajsf1_TrCer_1000_02_ddt_tracer.mat'],'tracklong')
%     trajs_conc = [trajs_conc tracklong];
% 2
%     load([folderin filesep 'trajsf2_TrCer_1000_02_ddt_tracer.mat'],'tracklong')
%     trajs_conc = [trajs_conc tracklong];
% 3
%      load([folderin filesep 'trajsf_TrCer_1000_04_ddt_tracer.mat'],'tracklong')
%      trajs_conc = [trajs_conc tracklong];

     load([folderin filesep 'trajsf_TrCer_1000_09_ddt_tracer.mat'],'tracklong')
     trajs_conc = [trajs_conc tracklong];
     2
     load([folderin filesep 'trajsf_TrCer_1000_10_ddt_tracer.mat'],'tracklong')
     trajs_conc = [trajs_conc tracklong];
     3

    Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);

end
trajs_conc = trajs_conc(Ine);
clear tracklong traj_fullg


%% Compute interpolant for mean flow 
dt = [4 6 8 10];
nbins = [20 21 22];
threshold = 1;
n=1;
power=1;

[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(trajs_conc,'Xf',dt,nbins,n,power,'x','cart');
[mYdt, mBydt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Yf',dt,nbins,n,power,'y','cart');
[mZdt, mBzdt, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Zf',dt,nbins,n,power,'z','cart');

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
[X,Y,Z]=ndgrid(linspace(-31,31,20),linspace(-45,36,21),linspace(-21,21,22));

Fx=griddedInterpolant(X,Y,Z,mXdts,'linear','none');
Fy=griddedInterpolant(X,Y,Z,mYdts,'linear','none');
Fz=griddedInterpolant(X,Y,Z,mZdts,'linear','none');

save([folderout filesep 'MeanFlow_Interpolant'],'Fx','Fy','Fz')

%% better way to plot center plane?

% Define grid for center plane
[x, y] = meshgrid(linspace(-20, 20, 30), linspace(-30, 30, 30));
z = -5.*ones(size(x)); % assuming the center plane is at z=0

% Evaluate interpolants on the grid
vx = Fx(x, y, z);
vy = Fy(x, y, z);
vz = Fz(x, y, z);

% Compute magnitude of velocity vector
v_mag = sqrt(vx.^2 + vy.^2 + vz.^2);

% Plot velocity vectors with colormap
figure;
quiver(x, -y, vx, vy, 1,'Color', 'k'); % plot velocity vectors
hold on;
surf(x,-y,v_mag, EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
colormap('jet');
c = colorbar;
c.Label.String = 'Velocity Magnitude';
caxis([min(v_mag(:)), max(v_mag(:))]); % set colorbar limits to the range of velocity magnitude
title('Velocity Vector Field');
xlabel('x');
ylabel('y');
axis equal;
hold off;
savefig_FC([folderout filesep 'meanfields_center' filesep name],8,6,'fig')
savefig_FC([folderout filesep 'meanfields_center' filesep name],8,6,'pdf')
stop
%ylim([-25 25])
%xlim([-15 15])

%% Tracer: compute the fluctuation of fluid velocity sigma_u = sqrt(<u'^2>)
dt = [4 6 8 10];
nbins = [20 21 22];
threshold = 1;
gridRange.x = [-30 30];
gridRange.y = [-45 35];
gridRange.z = [-20 20];

[gridsV,meanFieldsV,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,1,1,gridRange,1);
 [gridsVrms,meanFieldsVrms,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,1,2,gridRange,1);
 [gridsA,meanFieldsA,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,2,1,gridRange,1);
 [gridsArms,meanFieldsArms,trajs_conc] = meanFields(trajs_conc,Fs,dt,nbins,threshold,2,2,gridRange,1);

sigmau1.x = sqrt(mean(vertcat(trajs_conc.sVx).^2,"omitnan"))
sigmau1.y = sqrt(mean(vertcat(trajs_conc.sVy).^2,"omitnan"))
sigmau1.z = sqrt(mean(vertcat(trajs_conc.sVz).^2,"omitnan"))
sigmau1.all = sqrt((sigmau1.x^2+sigmau1.y^2+sigmau1.z^2)/3) % unit: mm/s

%%  Plot Mean Velocity CENTER

plane = 11;% CENTER
%plane = 2;% OFF-CENTER

mkdir([folderout filesep 'meanfields_center'])

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

savefig_FC([folderout filesep 'meanfields_center' filesep name],8,6,'fig')
savefig_FC([folderout filesep 'meanfields_center' filesep name],8,6,'pdf')

%%% Plot center plane at constant x -- couldn't make it work
% name = 'vertical_plane_cst_x';
% figure;
% 
% vel = reshape(normvel(:,:,50),[100,200]);
% Ys = reshape(gridsV.YY(:,:,50),[100,200]);
% Zs = reshape(gridsV.ZZ(50,:,:),[200,100]);
% surf(Ys,Ys,real(vel), EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
% %q=quiver(gridsV.ZZ(:,50,:),gridsV.YY(:,50,:),meanFieldsV.z(:,50,:),meanFieldsV.y(:,50,:),'r','LineWidth',1,'autoscalefactor',1);
% 
% %title(name,'Interpreter','latex')
% colorbar ; colormap jet; clim([0 200])
% axis equal
% % xlim([-13,10])
% % ylim([-30,20])
% xlabel('z(mm)')
% ylabel('y(mm)')
% view(0,90)
% stop
% % savefig_FC(['meanfields' filesep name],8,6,'fig')
% % savefig_FC(['meanfields' filesep name],8,6,'pdf')

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

