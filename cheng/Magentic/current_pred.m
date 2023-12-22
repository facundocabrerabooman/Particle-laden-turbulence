clc;clear;
close all

Colis_Param = {};

% Z6-Z1
Colis_Param.R = [16.1000   21.8000   15.4000   15.6000   22.1000   16.3000]*1e-2; % m
Colis_Param.N = [969   101   452   450   103   965];

% This Week (23/05/20223)
% wider region without sign change but the fluctuation is a bit high
Colis_Param.I =  [-2.43  -0.25   -0.15     1.79   0.25    2.0]*(0.5)
Colis_Param.Z = [-26.5000  -24.5000  -12.0000   12.0000   24.5000   26.5000]*1e-2; % m

grad = 1;
Z_TankCenter = 90;
delta_R= 0.01; % m

filelist = {'calibration_14_52_53';};
for k = 1:length(filelist)
    CalibData(:,:,k) = load(['D:\1-MageticData\', char(filelist(k)),'.txt']);
    z_Exp(:,k) = (CalibData(:,1,k)-Z_TankCenter)*1e-3; % m
end

color1 = colormap('jet');
close
color0 = color1(1:floor(size(color1,1)/length(filelist)):end,:);


%%
[B_pred_indv,B_pred0] = B_pred(Colis_Param,z_Exp(:,1));
Gz_pred = diff(B_pred0)./diff(z_Exp(:,1));
kk = find(abs(B_pred0) == min(abs(B_pred0)));
kk00 = find(abs(z_Exp) == min(abs(z_Exp)));
kk01 = find(abs(z_Exp-0.1) == min(abs(z_Exp-0.1)));
kk005= find(abs(z_Exp-0.05) == min(abs(z_Exp-0.05)));
Z_height = z_Exp(1)-z_Exp(kk);
% figure 
subplot(2,1,1)
plot(z_Exp(:,1),B_pred0,'LineWidth',2)
hold on 
axis tight
plot(xlim,[0 0],'r--')
plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
plot([z_Exp(kk00) z_Exp(kk00)],ylim,'b--')
plot([z_Exp(kk01) z_Exp(kk01)],ylim,'b--')
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$B/G$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
legend('Theory')
    
%%
subplot(2,1,2)
plot(z_Exp(1:end-1,1),Gz_pred,'LineWidth',2)
hold on
axis tight
plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
plot([z_Exp(kk00) z_Exp(kk00)],ylim,'b--')
plot([z_Exp(kk01) z_Exp(kk01)],ylim,'b--')
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$$\nabla_z B(G/m) $$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
legend('Theory')

%% Exp data compute
r = [0:0.005:0.05]'; % m 

if grad == 1     
    [X1,Y1] = meshgrid(z_Exp(:,1),[0:0.005:0.05]);
else 
    [X1,Y1] = meshgrid(z_Exp(1:end-1,1),[0:0.005:0.045]);
end
X1 = X1';
Y1 = Y1';
    
%% gravity compensation
% % 200micron cube
% B_rem = 1.24; % Tesla
% mass = 0.001e-3; % kg
% V = (0.4e-3)^3; % m^3
% miu0 = 4*pi*10^-7;% constant
% M = B_rem*V/miu0;

% 1mm sphere
B_rem = 1.17; % Tesla
mass = 0.01e-3; % kg
d = 1e-3; % m
V = 4/3*pi*(d/2)^3; % m^3
miu0 = 4*pi*10^-7; % constant
M = B_rem*V/miu0;

g = 9.8;
% g_eff = g-M*dBdz/10^4*1/mass;
Gz_pred0 = zeros(99,11);
for i = 1:11
    Gz_pred0(:,i) = Gz_pred;
end

g_eff = g-M*Gz_pred0/10^4*1/mass;

g_ratio = g_eff/g;

figure
contourf(Y1(1:end-1,:),X1(1:end-1,:),g_ratio,'ShowText','on',LevelStep=0.01)
hold on
plot(xlim,[0 0],'r--')
plot(xlim,[0.1 0.1],'r--')
colorbar
colormap('summer')
ylabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
xlabel('$r/m $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
title('$g_{eff}/g$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
