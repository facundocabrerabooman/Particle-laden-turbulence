function [ result ] = checkRotFilterDiscrete( ctraj, omegaFil, alphaFil, doAlpha )
%[ result ] = checkRotFilter( ctraj, omegaFil, diffOmegaFil, doMatrix )
%   Detailed explanation goes here

xyz={'x','y','z'};
omegaData =zeros(length(omegaFil),9);
alphaData=zeros(length(alphaFil),9);


which='PdLUR';
Props.technique='discrete';



%% velocity
for i=1:length(omegaFil)
  kernel.xw=1;
  kernel.vw=omegaFil(i);
  kernel.aw=2;
  
  Props.kernel=kernel;
  ctraj=VHT_angular(ctraj,'LUR',Props);
  
  for v=1:3
    vel=[ctraj.(['w' which xyz{v}])];%
    omegaData(i,0+v)=nanmean(vel);
    omegaData(i,3+v)=nanstd(vel);
    omegaData(i,6+v)=kurtosis(vel);
  end
end

if doAlpha
  kernel.xw=1;
  kernel.vw=omegaFil(1);
  kernel.aw=2;
  
  for i=1:length(alphaFil)
    
    kernel.aw=alphaFil(i);
    Props.kernel=kernel;
    
    ctraj=VHT_angular(ctraj,'LUR',Props);
    
    for v=1:3
      acc=[ctraj.(['a' which xyz{v}])];
      alphaData(i,0+v)=nanmean(acc);
      alphaData(i,3+v)=nanstd(acc);
      alphaData(i,6+v)=kurtosis(acc);
    end
  end
  
  
end







%% plotting
if ~doAlpha
  figure(5487837)
  
  subplot(3,1,1)
  hold off
  plot(omegaFil,omegaData(:,1:3),'x')
  ylabel('mean ang vel [rad/s]')
  legend({'<\omega_x>','<\omega_y>','<\omega_z>'},'Location','best')
  title('Discrete Technique, Vel, Gaussian kernels, btw: f_c=f_s/(1.42 W)')
  set(gca,'XMinorTick','on','YMinorTick','on')
  grid on
  
  subplot(3,1,2)
  hold off
  plot(omegaFil,omegaData(:,4:6),'x')
  ylabel('std(ang vel) [rad/s]')
  legend({'std(\omega_x)','std(\omega_y)','std(\omega_z)'},'Location','best')
  set(gca,'XMinorTick','on','YMinorTick','on')
  grid on
  
  subplot(3,1,3)
  hold off
  plot(omegaFil,omegaData(:,7:9),'x')
  xlabel('Filter Width [frames]')
  ylabel('Flatness(ang vel)')
  legend({'F(\omega_x)','F(\omega_y)','F(\omega_z)'},'Location','best')
  set(gca,'XMinorTick','on','YMinorTick','on')
  grid on
  
  
else
  figure(2784)
  subplot(3,1,1)
  hold off
  plot(alphaFil,alphaData(:,1:3),'x')
  ylabel('mean ang acc [rad/s^2]')
  legend({'<\alpha_x>','<\alpha_y>','<\alpha_z>'},'Location','best')
  title('Discrete Technique, Acc, Gaussian kernels, btw: f_c=f_s/(1.42 W)')
  set(gca,'XMinorTick','on','YMinorTick','on')
  grid on
  
  
  
  subplot(3,1,2)
  hold off
  plot(alphaFil,alphaData(:,4:6),'x')
  ylabel('std(ang acc) [rad/s^2]')
  legend({'std(\alpha_x)','std(\alpha_y)','std(\alpha_z)'},'Location','best')
  set(gca,'XMinorTick','on','YMinorTick','on')
  grid on
  
  
  subplot(3,1,3)
  hold off
  semilogy(alphaFil,alphaData(:,7:9),'x')
  xlabel('Filter Width [frames]')
  ylabel('Flatness(ang acc)')
  legend({'F(\alpha_x)','F(\alpha_y)','F(\alpha_z)'},'Location','best')
  set(gca,'XMinorTick','on','YMinorTick','on')
  grid on
end

result.acc=alphaData;
result.filterAcc=alphaFil;
result.vel=omegaData;
result.filterVel=omegaFil;

end

