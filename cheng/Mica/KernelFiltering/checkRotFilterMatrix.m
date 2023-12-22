function [ result ] = checkRotFilterMatrix( ctraj, omegaFil, diffOmegaFil )
%[ result ] = checkRotFilterMatrix( ctraj, omegaFil, diffOmegaFil)
%   Detailed explanation goes here

xyz={'x','y','z'};
omegaData =zeros(length(omegaFil),9);
accOmgData=zeros(length(diffOmegaFil),9);


which='wPmLRR';
Props.doRunningAvg=false; % no running average
Props.maxAvg=10;
Props.technique='matrix'; % 'discrete'
ctraj=VHT_angular(ctraj,'LRR',Props);

KernelParam.doPos=true;
KernelParam.doVel=false;
KernelParam.doAcc=false;

KernelParam.prePos='k';
KernelParam.preVel='kd';
%% smoothing
for i=1:length(omegaFil)
  KernelParam.xw=omegaFil(i);
  
  ctraj=VHT_Kernel(ctraj,which,KernelParam);
  
  for v=1:3
    vel=[ctraj.(['k' which xyz{v}])];%
    omegaData(i,0+v)=nanmean(vel);
    omegaData(i,3+v)=nanstd(vel);
    omegaData(i,6+v)=kurtosis(vel);
  end
end

%% deriving
KernelParam.doPos=false;
KernelParam.doVel=true;

for i=1:length(diffOmegaFil)
  
  KernelParam.vw=diffOmegaFil(i);
  
  ctraj=VHT_Kernel(ctraj,which,KernelParam);
  
  for v=1:3
    acc=[ctraj.(['kd' which xyz{v}])];
    accOmgData(i,0+v)=nanmean(acc);
    accOmgData(i,3+v)=nanstd(acc);
    accOmgData(i,6+v)=kurtosis(acc);
  end
end




%% plotting


figure(587897)
title('Matrix Technique, Gaussian kernels, btw: f_c=f_s/(1.42 W)')
subplot(3,2,1)
hold off
plot(omegaFil,omegaData(:,1:3),'x')
% xlabel('Filter Width [frames]')
ylabel('mean ang vel [rad/s]')
legend({'<\omega_x>','<\omega_y>','<\omega_z>'},'Location','best')
title('Matrix Technique, Gaussian kernels, btw: f_c=f_s/(1.42 W)')
grid on 
set(gca,'XMinorTick','on','YMinorTick','on')

subplot(3,2,2)
hold off
plot(diffOmegaFil,accOmgData(:,1:3),'x')
% xlabel('Filter Width [frames]')
ylabel('mean ang acc [rad/s^2]')
legend({'<\alpha_x>','<\alpha_y>','<\alpha_z>'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on')
grid on 

subplot(3,2,3)
hold off
plot(omegaFil,omegaData(:,4:6),'x')
% xlabel('Filter Width [frames]')
ylabel('std(ang vel) [rad/s]')
legend({'std(\omega_x)','std(\omega_y)','std(\omega_z)'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on')
grid on 

subplot(3,2,4)
hold off
plot(diffOmegaFil,accOmgData(:,4:6),'x')
% xlabel('Filter Width [frames]')
ylabel('std(ang acc) [rad/s^2]')
legend({'std(\alpha_x)','std(\alpha_y)','std(\alpha_z)'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on')
grid on 

subplot(3,2,5)
hold off
plot(omegaFil,omegaData(:,7:9),'x')
xlabel('Filter Width [frames]')
ylabel('Flatness(ang vel)')
legend({'F(\omega_x)','F(\omega_y)','F(\omega_z)'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on')
grid on 

subplot(3,2,6)
hold off
semilogy(diffOmegaFil,accOmgData(:,7:9),'x')
xlabel('Filter Width [frames]')
ylabel('Flatness(ang acc)')
legend({'F(\alpha_x)','F(\alpha_y)','F(\alpha_z)'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on')
grid on 

result.acc=accOmgData;
result.filterAcc=diffOmegaFil;
result.vel=omegaData;
result.filterVel=omegaFil;

end

