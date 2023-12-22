clc;clear;
close all

Colis_Param = {};

% Z6-Z1
Colis_Param.R = [16.3000   22.1000   15.6000   15.4000   21.8000   16.1000]*1e-2; % m
Colis_Param.N = [965   103   450   452   101   969];

% 23/05/2023
% wider region without sign change but the fluctuation is a bit high
Colis_Param.I = [-2.43 -0.25 -0.15 1.79 0.25 2.0]*2; %calibration value
Colis_Param.Z = [-26.5000  -24.5000  -12.0000   12.0000   24.5000   26.5000]*1e-2; % m

CalibData = load('D:\1-MageticData\calibration_14_52_53.txt');
% FacusGz = xlsread('D:\Cheng\Data\Calib_Magnetic\Facus Dataset_ adapted from figure in thesis\FacusDataset-Gradient.xlsx');
% FacusB = xlsread('D:\Cheng\Data\Calib_Magnetic\Facus Dataset_ adapted from figure in thesis\FacusDataset-Magnetic.xlsx');

Z_TankCenter = 90;
z_Exp = (CalibData(:,1)-Z_TankCenter)*1e-3; % m
Bz_Exp = CalibData(:,4); % Gauss
% z_Exp = (-0.2:0.001:0.2)'; % m

%% 
% [B_pred_indv,B_pred0] = B_pred(Colis_Param,z_Exp);
% Gz_pred = diff(B_pred0)./diff(z_Exp);
% kk = find(abs(B_pred0) == min(abs(B_pred0)));
% kk00 = find(abs(z_Exp) == min(abs(z_Exp)));
% kk01 = find(abs(z_Exp-0.1) == min(abs(z_Exp-0.1)));
% Z_height = z_Exp(1)-z_Exp(kk)
% %     figure 
% subplot(2,1,1)
% plot(z_Exp,B_pred0)
% axis tight
% hold on
% plot(xlim,[0 0],'r--')
% plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
% plot([z_Exp(kk00) z_Exp(kk00)],ylim,'b--')
% plot([z_Exp(kk01) z_Exp(kk01)],ylim,'b--')
% xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
% ylabel('$B/G$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
% subplot(2,1,2)
% plot(z_Exp(1:end-1),Gz_pred)
% hold on 
% plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
% plot([z_Exp(kk00) z_Exp(kk00)],ylim,'b--')
% plot([z_Exp(kk01) z_Exp(kk01)],ylim,'b--')
% axis tight
% xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
% ylabel('$$\Delta B(G/m) $$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
% % hold off
% 
% 
% diff_Gz = max(Gz_pred(kk01:kk00))-min(Gz_pred(kk01:kk00))
% mean_Gz = mean(Gz_pred(kk01:kk00))
% fluc = abs(diff_Gz/mean_Gz)
% % disp([diff_Gz,mean_Gz,fluc])

%% compare exp and theory prediction
[B_pred_indv,B_pred0] = B_pred(Colis_Param,z_Exp);
Gz_pred = diff(B_pred0)./diff(z_Exp);
kk = find(abs(B_pred0) == min(abs(B_pred0)));
% kk00 = find(abs(z_Exp) == min(abs(z_Exp)));
% kk01 = find(abs(z_Exp-0.1) == min(abs(z_Exp-0.1)));
kk0015 = find(abs(z_Exp+0.015) == min(abs(z_Exp+0.015)));
kk0075 = find(abs(z_Exp-0.075) == min(abs(z_Exp-0.075)));

Z_height = z_Exp(1)-z_Exp(kk)
figure 
subplot(2,1,1)
plot(z_Exp,B_pred0,'LineWidth',2)
hold on 
plot(z_Exp,Bz_Exp,'ro')
axis tight
plot(xlim,[0 0],'r--')
plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
plot([z_Exp(kk0015) z_Exp(kk0015)],ylim,'g--')
plot([z_Exp(kk0075) z_Exp(kk0075)],ylim,'b--')
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$B/G$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
legend('Theory','EXP')
subplot(2,1,2)
plot(z_Exp,Bz_Exp-B_pred0,'bo')
axis tight
hold on
plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
plot([z_Exp(kk0015) z_Exp(kk0015)],ylim,'g--')
plot([z_Exp(kk0075) z_Exp(kk0075)],ylim,'b--')
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$$\Delta (G) $$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
%%
Gz_Exp = diff(Bz_Exp)./diff(z_Exp);
figure 
subplot(2,1,1)
plot(z_Exp(1:end-1),Gz_pred,'LineWidth',2)
hold on
plot(z_Exp(1:end-1),Gz_Exp,'ro')
axis tight
plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
plot([z_Exp(kk0015) z_Exp(kk0015)],ylim,'g--')
plot([z_Exp(kk0075) z_Exp(kk0075)],ylim,'b--')
legend('Theory','EXP')
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$$\nabla_z B(G/m) $$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
subplot(2,1,2)
plot(z_Exp(1:end-1),Gz_Exp-Gz_pred,'bo')
hold on 
axis tight
plot([z_Exp(kk) z_Exp(kk)],ylim,'r--')
plot([z_Exp(kk0015) z_Exp(kk0015)],ylim,'g--')
plot([z_Exp(kk0075) z_Exp(kk0075)],ylim,'b--')
xlabel('$z/m$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$$\Delta(G/m)$$',FontSize=20,FontName='Times New Roman',Interpreter='latex')

disp(['fluctuation = ' num2str(max(Gz_Exp(kk0075:kk0015))-min(Gz_Exp(kk0075:kk0015)))])
disp(['Mean Gradient = ' num2str(mean(Gz_Exp(kk0075:kk0015)))])
disp(['Fluc/Mean = ' num2str((max(Gz_Exp(kk0075:kk0015))-min(Gz_Exp(kk0075:kk0015)))/mean(Gz_Exp(kk0075:kk0015)))])

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

deltaG = M*mean(Gz_pred)/10^4*1/mass/g

a.Z = fliplr(Colis_Param.Z);
a.I = fliplr(Colis_Param.I)
1-deltaG



