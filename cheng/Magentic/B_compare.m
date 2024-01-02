clc;clear;
close all
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));
mycolormap1 = mycolor('#0000B2','#FFFFFF','#B10000');

%% B field case
% 1: micro-gravity
% 2: hyper-gravity
% 3: constant B (normal gravity)
ifconfig = 3; 

%% position values in the stage coordinate
Z_ColisCenter =90; %center of the coli sets instead of the geometric center of the tank

image_bot = 95;
image_top = 171;

cam_pos = 133; % find with set_0.PY, equals to the center of the images

traj_bot = 95; % roughly equals to the postion of the bottom of the image
traj_top = 167; % roughly equals to the postion of the end of the tube
% length of the trajectories is around 68~70 mm 

%% position values in the colis' coordinate
image_bot = (image_bot - Z_ColisCenter)*1e-3;
image_top = (image_top - Z_ColisCenter)*1e-3;

cam_pos = (cam_pos - Z_ColisCenter)*1e-3;

traj_bot = (traj_bot - Z_ColisCenter)*1e-3;
traj_top = (traj_top - Z_ColisCenter)*1e-3;

%%

Colis_Param = {};

% Z6-Z1
Colis_Param.R = [16.3000   22.1000   15.6000   15.4000   21.8000   16.1000]*1e-2; % m
Colis_Param.N = [965   103   450   452   101   969];

%% micro gravity (23/05/20223)
% wider region without sign change but the fluctuation is a bit high
if ifconfig == 1
    Colis_Param.I =  [2.0  0.25  1.79  -0.15  -0.25  -2.43]*2;
    Colis_Param.Z = [26.5000   24.5000   12.0000  -12.0000  -24.5000  -26.5000]*1e-2; % m
    
    grad = 1;
    delta_R= 0.005; % m
    
    % micro gravity calibration
    fpath = 'D:\1-MageticData-Microgravity\';
    filelist = {'calibration_14_52_53';...
                'calibration_09_44_18';...
                'calibration_15_18_20';...
                'calibration_10_15_54';...
                'calibration_10_31_56';...
                'calibration_10_47_06';...
                'calibration_11_02_08';...
                'calibration_11_17_35';...
                'calibration_11_33_10';...
                'calibration_11_47_59';
                'calibration_13_54_18'};
    
    LevelStep_F = 0.05;
    LevelStep_g = 0.02;
end
%% Hyper gravity (05/10/20223)
if ifconfig == 2
    Colis_Param.I =  [-0.85,  0.48,  0.25,   0.50,   0.75,   0.608];
    Colis_Param.Z = [26.5000   24.5000   12.0000  -12.0000  -24.5000  -26.5000]*1e-2; % m
    
    grad = 1;
    delta_R= 0.005; % m
    
    % micro gravity calibration
    fpath = 'D:\2-MageticData-HyperGravity\';
    filelist = {'calibration_16_44_18';...
                'calibration_17_35_45';...
                'calibration_17_57_08';...
                'calibration_18_14_13';...
                'calibration_18_30_26';...
                'calibration_18_47_25';...
                'calibration_19_04_09';...
                'calibration_19_23_16';...
                'calibration_19_40_39';...
                'calibration_09_45_58';...
                'calibration_10_02_02'};

    LevelStep_F = 0.05;
    LevelStep_g = 0.02;

end

%% constant B normal gravity
if ifconfig == 3
    Colis_Param.I =  [-0.10,  0.75,  0.06,   0.22,   0.75,   -0.39];
    Colis_Param.Z = [26.5000   24.5000   12.0000  -12.0000  -24.5000  -26.5000]*1e-2; % m
    
    grad = 1;
    delta_R= 0.005; % m
    
    % micro gravity calibration
    fpath = 'D:\3-MageticData-ConstantB\';
    filelist = {'calibration_10_39_51';...
                'calibration_10_54_56';...
                'calibration_11_09_38';...
                'calibration_11_24_21';...
                'calibration_11_38_49';...
                'calibration_11_53_10';...
                'calibration_12_07_28';...
                'calibration_13_03_57';...
                'calibration_13_20_01';...
                'calibration_13_34_13';...
                'calibration_13_48_33'};
    
    LevelStep_F = 20;
    LevelStep_g = 0.005;

end

%%
for k = 1:length(filelist)
    CalibData(:,:,k) = load([fpath, char(filelist(k)),'.txt']);
    z_Exp(:,k) = (CalibData(:,1,k)-Z_ColisCenter)*1e-3; % m
    Bx_Exp(:,k)  = CalibData(:,2,k); % Gauss
    By_Exp(:,k)  = CalibData(:,3,k); % Gauss
    Bz_Exp(:,k)  = CalibData(:,4,k); % Gauss

    Br_Exp(:,k)  = sqrt(Bx_Exp(:,k).^2 + By_Exp(:,k).^2); % Gauss

    B_EXP(:,k)  = sqrt(Bx_Exp(:,k).^2 + By_Exp(:,k).^2 + Bz_Exp(:,k).^2);
end
delta_z = z_Exp(2,1)-z_Exp(1,1);

color0 = mycolormap1(1:floor(size(mycolormap1,1)/length(filelist)):end,:);
%%
[B_pred_indv,B_pred0] = B_pred(Colis_Param,z_Exp(:,1));
Gz_pred = diff(B_pred0)./diff(z_Exp(:,1));
kk = find(abs(B_pred0) == min(abs(B_pred0)));
kk_cam_pos = find(abs(z_Exp-cam_pos) == min(abs(z_Exp-cam_pos)));
Z_height = z_Exp(1)-z_Exp(kk);
figure 
subplot(2,1,1)
plot(z_Exp(:,1),B_pred0,'LineWidth',2)
hold on 
for k = 1:length(filelist)
    plot(z_Exp(:,k),Bz_Exp(:,k),'o','LineWidth',2,Color=color0(k,:))
end
axis tight
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$B/G$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
legend('Theory','r=0',...
    'r=+0.005','r=+0.01', ...
    'r=+0.015','r=+0.02', ...
    'r=+0.025','r=+0.03', ...
    'r=+0.035','r=+0.04', ...
    'r=+0.045','r=+0.05','AutoUpdate','off')
yline(0,'r--')
xline(z_Exp(kk),'r--')
xline(cam_pos,'k--')

patch([image_bot image_bot image_top image_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#FFFFFF'),EdgeColor = 'b')
alpha(0)
patch([traj_bot traj_bot traj_top traj_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#a155b9'),EdgeColor = 'none')
alpha(0.2)


if grad == 1
    [dBzdr_Exp, dBzdz_Exp] = gradient(Bz_Exp,delta_R,delta_z);
    [dBrdr_Exp, dBrdz_Exp] = gradient(Br_Exp,delta_R,delta_z);

else   
    dBzdz_Exp = diff(Bz_Exp,1,1)./diff(z_Exp(:,1));
    dBzdr_Exp = diff(Bz_Exp,1,2)./delta_R;
    dBrdz_Exp = diff(Br_Exp,1,1)./diff(z_Exp(:,1));
    dBrdr_Exp = diff(Br_Exp,1,2)./delta_R;
end
    
%%
subplot(2,1,2)
plot(z_Exp(1:end-1,1),Gz_pred,'LineWidth',2)
hold on
for k = 1:length(filelist)
    if grad == 1
        plot(z_Exp(:,k),dBzdz_Exp(:,k),'o','LineWidth',2,Color=color0(k,:))
    else
        plot(z_Exp(1:end-1,k),dBzdz_Exp(:,k),'o','LineWidth',2,Color=color0(k,:))
    end
end
axis tight
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$$\nabla_z B(G/m) $$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
legend('Theory','r=0',...
    'r=+0.005','r=+0.01', ...
    'r=+0.015','r=+0.02', ...
    'r=+0.025','r=+0.03', ...
    'r=+0.035','r=+0.04', ...
    'r=+0.045','r=+0.05','AutoUpdate','off')

xline(z_Exp(kk),'r--')
xline(cam_pos,'k--')

patch([image_bot image_bot image_top image_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#FFFFFF'),EdgeColor = 'b')
alpha(0)
patch([traj_bot traj_bot traj_top traj_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#a155b9'),EdgeColor = 'none')
alpha(0.2)


%% Exp data compute
r = [0:delta_R:delta_R*(length(filelist)-1)]'; % m 

if grad == 1  
    dBdz = (Br_Exp.*dBrdz_Exp + Bz_Exp.*dBzdz_Exp)./B_EXP;
    dBdr = (Br_Exp.*dBrdr_Exp + Bz_Exp.*dBzdr_Exp)./B_EXP;
else
    dBdz = (Br_Exp(1:end-1,:).*dBrdz_Exp + Bz_Exp(1:end-1,:).*dBzdz_Exp)./B_EXP(1:end-1,:);
    dBdr = (Br_Exp(:,1:end-1).*dBrdr_Exp + Bz_Exp(:,1:end-1).*dBzdr_Exp)./B_EXP(:,1:end-1);
end

    ratio_zz = dBdz/dBdz(kk_cam_pos(1),1);
%     ratio_rz = 1- dBdr/dBdz(kk003(1),1);
    ratio_rz = dBdr/dBdz(kk_cam_pos(1),1);

if grad == 1     
    [X1,Y1] = meshgrid(z_Exp(:,1),[0:delta_R:delta_R*(length(filelist)-1)]);
else 
    [X1,Y1] = meshgrid(z_Exp(1:end-1,1),[0:0.005:0.045]);
end
X1 = X1';
Y1 = Y1';
figure
subplot(1,2,1)
% pcolor(Y1,X1,ratio_zz);shading interp
contourf(Y1,X1,ratio_zz,'ShowText','on',LevelStep=LevelStep_F)
hold on
yline(image_top,'b--')
yline(image_bot,'b--')
yline(traj_top,'r--')
yline(traj_bot,'r--')
colorbar
ylabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
xlabel('$r/m $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
title('$F_z(r,z)/F_z(0,z_{cam})$',FontSize=20,FontName='Times New Roman',Interpreter='latex')

subplot(1,2,2)
%     pcolor(Y1,X1,ratio_rz);shading interp
contourf(Y1,X1,ratio_rz,'ShowText','on',LevelStep=LevelStep_F)
hold on
yline(image_top,'b--')
yline(image_bot,'b--')
yline(traj_top,'r--')
yline(traj_bot,'r--')
colormap("summer")
colorbar
ylabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
xlabel('$r/m $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
title('$F_r(r,z)/F_z(0,z_{cam})$',FontSize=20,FontName='Times New Roman',Interpreter='latex')

% ratio_rz_fix = dBdr./dBdz;
% subplot(1,3,3)
% %     pcolor(Y1,X1,ratio_rz_fix);shading interp
% contourf(Y1,X1,ratio_rz_fix,'ShowText','on',LevelStep=0.02)
% hold on
% plot(xlim,[0 0],'r--')
% plot(xlim,[0.1 0.1],'r--')
% colorbar
% caxis([-1 1])
% ylabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
% xlabel('$r/m $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
% title('$F_r(r,z)/F_z(r,z)$',FontSize=20,FontName='Times New Roman',Interpreter='latex')

    
%% gravity compensation
% 400micron cube
% B_rem = 1.24; % Tesla
% mass = 0.001e-3; % kg
% V = (0.4e-3)^3; % m^3
% miu0 = 4*pi*10^-7;% constant
% M = B_rem*V/miu0;

% 1mm sphere
B_rem = 1.17; % Tesla
mass = 0.0048e-3; % kg
d = 1e-3; % m
V = 4/3*pi*(d/2)^3; % m^3
miu0 = 4*pi*10^-7; % constant
M = B_rem*V/miu0;

g = 9.8;
g_eff = g-M*dBdz/10^4*1/mass;
g_ratio = g_eff/g; 

figure
contourf(Y1,X1,g_ratio,'ShowText','on',LevelStep=LevelStep_g)
hold on
yline(image_top,'b--')
yline(image_bot,'b--')
yline(traj_top,'r--')
yline(traj_bot,'r--')
colorbar
colormap("summer")
ylabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
xlabel('$r/m $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
title('$g_{eff}/g$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
