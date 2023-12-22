%% Definitions
% Facu's lsq fit
clear all


%to change direction of grad B; change sign of A!
A = 0/1e4; %gradB
%B = 5/1e4; % B(0)
B = 1.25/1e4; % B(0)

%%% data from fit
Colis_Param.R = [16.1000   21.8000   15.4000   15.6000   22.1000   16.3000]*1e-2; % m
Colis_Param.N = [969   101   452   450   103   965];
Colis_Param.Z = [-26.5000  -24.5000  -12.0000   12.0000   24.5000   26.5000]*1e-2;

R1 = Colis_Param.R(1);d1 = Colis_Param.Z(1);N1 = Colis_Param.N(1);
R2 = Colis_Param.R(2);d2 = Colis_Param.Z(2);N2 = Colis_Param.N(2);
R3 = Colis_Param.R(3);d3 = Colis_Param.Z(3);N3 = Colis_Param.N(3);
R4 = Colis_Param.R(4);d4 = Colis_Param.Z(4);N4 = Colis_Param.N(4);
R5 = Colis_Param.R(5);d5 = Colis_Param.Z(5);N5 = Colis_Param.N(5);
R6 = Colis_Param.R(6);d6 = Colis_Param.Z(6);N6 = Colis_Param.N(6);


%%% equation 
% 1a1b = d1; 1c1d=d2; 2a2b = d3...
func = @(x,z) 2*pi*10^(-7)*( R1.^2*N1*x(1)./((z-d1).^2 + R1.^2).^(1.5) + N2*R2.^2*x(2)./((z-d2).^2 + R2.^2).^(1.5) + R3.^2*N3*x(3)./((z-d3).^2 + R3.^2).^(1.5) + R4.^2*N4*x(4)./((z-d4).^2 + R4.^2).^(1.5) + R5.^2*N5*x(5)./((z-d5).^2 + R5.^2).^(1.5)+ R6.^2*N6*x(6)./((z-d6).^2 + R6.^2).^(1.5)) ;

options = optimoptions('lsqcurvefit');
options.OptimalityTolerance = 1e-25;
options.MaxIterations = 1e20;
options.MaxFunctionEvaluations = 1e25;
options.FunctionTolerance = 1e-25;
options.StepTolerance = 1e-22;

%toped = d1;
%topeu = d2;
topeu = -0.05;
toped = 0.1;
z = (topeu:0.0001:toped);
z_plot = (-0.20:0.0001:0.20);
y = A*z + B;
y_plot = A*z_plot + B;
IC =  [1 1 1 1 1 1];
Lb = -[3 3 3 3 3 3];
Ub =  [3 3 3 3 3 3];
%exitflag,output
[x,resnorm] = lsqcurvefit(func,IC,z,y,Lb,Ub,options);
%%
figure %clf , 
a=rand(1,3); % color chosen randomly 
subplot(1,2,1);plot(z_plot,y_plot*10000,'.r'),grid on, hold on, plot(z_plot,func(x,z_plot)*10000,'-b')
xline(topeu,'color',a,'Linewidth',1);xline(toped,'color',a,'Linewidth',1);
yline(0);xline(0);
yline(5);
yline(-5);
ylabel('B [G]')
xlabel('dist to center [m]')
title('B')
subplot(1,2,2);plot(z_plot(1:end-1),diff(func(x,z_plot)*10000)./1e-4,'color',a,'Linewidth',3), xlim([-0.2 0.2]); 
hold on
xline(topeu,'color',a,'Linewidth',1);xline(toped,'color',a,'Linewidth',1);
%ylim([A*1.02e4 A*0.98e4])
% ylim([200 300])
ylabel('\nablaB [G/m]')
xlabel('dist to center [m]')
title('gradB')

gradB = diff(func(x,z)*10000)./1e-4;
%Jdown = find(ismembertol(z,-tope,0.001)==1);
%Jup = find(ismembertol(z,tope,0.001)==1);
mean(gradB)
std(gradB)

std(gradB)/mean(gradB)
