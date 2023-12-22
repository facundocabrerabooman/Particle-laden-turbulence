clc;clear;
close all

%%
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

% micro gravity
% Colis_Param.I =  -[2.0  0.25  1.79  -0.15  -0.25  -2.43]*2; % calibration value
% reverse the current doesn't change the sign of the magentic force, the force
% direction is only related to the magnitude of the Magenetic field(i.e: abs(B))

% hyper gravity
% Colis_Param.I = [-0.85,  0.48,  0.25,   0.50,   0.75,   0.608]; % calibration value, 1.1g

% constant B field
% calibrated at B(0) = 1.25 G constant, ~ 1.0037g, currents = [-0.3968    3.0000    0.2472    0.8640    3.0000   -1.5416]*0.25
% experiment performed at B(0)= 5 G, ~1.015g,    currents = [-0.3968    3.0000    0.2472    0.8640    3.0000   -1.5416]
Colis_Param.I = [-0.3968    3.0000    0.2472    0.8640    3.0000   -1.5416]*0.25; 

% Colis_Param.I =  [-1.7 0.96 0.51 0.99 1.5 1.215];
Colis_Param.Z = [26.5000   24.5000   12.0000  -12.0000  -24.5000  -26.5000]*1e-2; % m


%% micro gravity (23/05/20223)
% CalibData = load('D:\1-MageticData-Microgravity\calibration_14_52_53.txt');
% 
% Z_ColisCenter = 90; %center of the coli sets instead of the geometric center of the tank
% z_Exp = (CalibData(:,1)-Z_ColisCenter)*1e-3; % m
% Bz_Exp = CalibData(:,4); % Gauss
% % z_Exp = (-0.2:0.001:0.2)'; % m

%% Hyper gravity (05/10/20223)

% CalibData = load('D:\2-MageticData-HyperGravity\calibration_16_44_18.txt');
% 
% Z_ColisCenter = 90; %center of the coli sets instead of the geometric center of the tank
% z_Exp = (CalibData(:,1)-Z_ColisCenter)*1e-3; % m
% Bz_Exp = CalibData(:,4); % Gauss
z_Exp = (-0.09:0.001:0.13)'; % m

%% compare exp and theory prediction
[B_pred_indv,B_pred0] = B_pred(Colis_Param,z_Exp);
Gz_pred0 = diff(B_pred0)./diff(z_Exp);
Gz_pred = abs(Gz_pred0);

kk = find(abs(B_pred0) == min(abs(B_pred0)));
subplot(2,1,1) 
% figure
plot(z_Exp,B_pred0,'LineWidth',2);hold on
axis tight
% yline(0,'r--')
xline(z_Exp(kk),'r--')
xline(cam_pos,'k--')

patch([image_bot image_bot image_top image_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#FFFFFF'),EdgeColor = 'b')
alpha(0)
patch([traj_bot traj_bot traj_top traj_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#a155b9'),EdgeColor = 'none')
alpha(0.2)
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$B/G$',FontSize=20,FontName='Times New Roman',Interpreter='latex')


%%
subplot(2,1,2) 
% figure
plot(z_Exp(1:end-1),Gz_pred0,'LineWidth',2)
hold on
% yline(0,'r--')
xline(z_Exp(kk),'r--')
xline(cam_pos,'k--')

patch([image_bot image_bot image_top image_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#FFFFFF'),EdgeColor = 'b')
alpha(0)
patch([traj_bot traj_bot traj_top traj_top], [min(ylim) max(ylim) max(ylim) min(ylim)],mycolor('#a155b9'),EdgeColor = 'none')
alpha(0.2)
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$$\nabla_z B(G/m) $$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
axis tight
%%
% figure
% colors = ['r','g','b','m','y','k'];
% for i = 1:size(B_pred_indv,2)
%     plot(z_Exp,B_pred_indv(:,i),'LineWidth',2, 'Color',colors(i));hold on
% end
% legend('CC0','CC1','CC2','CC3','CC4','CC5')


%%
k = find(z_Exp>traj_bot & z_Exp<traj_top);
num1 = max(Gz_pred(k));
num2 = min(Gz_pred(k));
num3 = (max(Gz_pred(k))-min(Gz_pred(k))) /mean(Gz_pred(k));

disp(['max B grad = ' num2str(num1)])
disp(['min B grad = ' num2str(num2)])
disp(['inhomogeneous = ' num2str(num3)])
%% gravity compensation
% 1mm sphere
B_rem = 1.17; % Tesla
mass = 0.0048e-3; % kg
d = 1e-3; % m
V = 4/3*pi*(d/2)^3; % m^3
miu0 = 4*pi*10^-7; % constant
M = B_rem*V/miu0;

g = 9.8;

deltaG = M*Gz_pred/10^4*1/mass/g;
mean(deltaG);
% 

Colis_Param.I

if abs(B_pred0(z_Exp==0.1))>abs(B_pred0(z_Exp==0))
    sign0 = -1;
elseif abs(B_pred0(z_Exp==0.1))<abs(B_pred0(z_Exp==0))
    sign0 = +1;
end
1+sign0*mean(deltaG)


