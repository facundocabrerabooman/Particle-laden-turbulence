clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Volumes/landau1/Tracers/ddt_moredata';
folderout = folderin;
cd(folderin)

Fs=2990; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Concatenate data

if pi==pi
    trajs_conc = [];

    load('trajs_TrCer_1000_11_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong

    load('trajs_TrCer_1000_13_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong

    load('trajs_TrCer_1000_14_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong

    load('trajs_TrCer_1000_15_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong

    load('trajsf_TrCer_1000_27_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong

    load('trajsf_TrCer_1000_28_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    
    load('trajsf_TrCer_1000_29_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    
    load('trajsf_TrCer_1000_30_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    
    load('trajsf_TrCer_1000_31_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    
    load('trajsf1_TrCer_1000_32_ddt_tracers.mat','tracklong1')
    trajs_conc = [trajs_conc tracklong1];
    clear tracklong1
    
    load('trajsf2_TrCer_1000_32_ddt_tracers.mat','tracklong2')
    trajs_conc = [trajs_conc tracklong2];
    clear tracklong2

    Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);

end
trajs_conc = trajs_conc(Ine);
clear tracklong traj_ddt

%% Tracer: compute the fluctuation of fluid velocity sigma_u = sqrt(<u'^2>)
dt = [4 6 8 10];
nbins = [20 21 22];
threshold = 10;
gridRange.x = [-40 40];
gridRange.y = [-40 40];
gridRange.z = [-40 40];

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
mkdir('meanfields')

normvel = sqrt((meanFieldsV.x).^2 + (meanFieldsV.y).^2 + (meanFieldsV.z).^2);

%%% Plot center plane at constant z
name = 'vertical_plane_cst_z';
figure;

vel = reshape(normvel(:,:,12),[25,50]);
Ys = reshape(gridsV.YY(:,:,12),[25,50]);
Xs = reshape(gridsV.XX(:,:,12),[25,50]);
surf(Ys,-Xs,real(vel), EdgeColor = 'none',FaceAlpha = 0.7,FaceColor='interp'); hold on
q=quiver(gridsV.YY(:,:,12),-gridsV.XX(:,:,12),meanFieldsV.y(:,:,12),-meanFieldsV.x(:,:,12),'k','LineWidth',1);

%title(name,'Interpreter','latex')
cb = colorbar ; colormap jet; clim([0 150]); 
axis equal
ylim([-30,30])
xlim([-15,15])
xlabel('x(mm)')
ylabel('y(mm)')
ylabel(cb,'$|\mathbf{v}|$ (mm/s)','Interpreter','latex','Rotation',270,'FontSize',20)
view(0,90)
box; grid off

savefig_FC(['meanfields' filesep name],8,6,'fig')
savefig_FC(['meanfields' filesep name],8,6,'pdf')

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

