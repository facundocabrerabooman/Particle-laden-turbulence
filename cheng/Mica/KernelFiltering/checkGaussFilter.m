function [ result ] = checkGaussFilter( ctraj, velFil, accFil )
%CHECKACCFILTER Summary of this function goes here
%   Detailed explanation goes here

velData=zeros(length(velFil),9);
accData=zeros(length(accFil),9);

KernelParam.doPos=false;
KernelParam.doVel=true;
KernelParam.doAcc=true;
KernelParam.reNan=  -1;% deactivate NaNs with a dilate


dt=ctraj(1).t(2)-ctraj(1).t(1);
xyz={'x','y','z'};

if isempty(intersect(fieldnames(ctraj),'interpolated'))
  intp=false(size([ctraj.('x')]));
else
  intp=[ctraj.interpolated];
end;

%%

for i=1:min(length(velFil),length(accFil))
  
  KernelParam.xw=1;
  KernelParam.xL=1;
  KernelParam.vw=velFil(i);
  KernelParam.aw=accFil(i);
  KernelParam.dt=dt;
  
  trjtemp=VHT_Kernel(ctraj,'',KernelParam);
  %   useAcc=~bwmorph(intp>0,'dilate',2);
  %   useVel=~bwmorph(intp>0,'dilate',3);
  
  for v=1:3
    vel=[trjtemp.(['kv' xyz{v}])];
    velData(i,0+v)=nanmean(vel);
    velData(i,3+v)=nanstd(vel);
    velData(i,6+v)=kurtosis(vel);
    
    acc=[trjtemp.(['ka' xyz{v}])];
    accData(i,0+v)=nanmean(acc);
    accData(i,3+v)=nanstd(acc);
    accData(i,6+v)=kurtosis(acc);
  end
end

for i=(min(length(velFil),length(accFil))+1):max(length(velFil),length(accFil))
  
  if length(velFil)>length(accFil)
    KernelParam.vw=velFil(i);
    KernelParam.doPos=false;
    KernelParam.doVel=true;
    KernelParam.doAcc=false;
  else
    KernelParam.aw=accFil(i);
    KernelParam.doPos=false;
    KernelParam.doVel=false;
    KernelParam.doAcc=true;
  end
  
  trjtemp=VHT_Kernel(ctraj,'',KernelParam);
  
  if length(velFil)>length(accFil)
    for v=1:3
      vel=[trjtemp.(['kv' xyz{v}])];
      velData(i,0+v)=nanmean(vel);
      velData(i,3+v)=nanstd(vel);
      velData(i,6+v)=kurtosis(vel);
    end
    
  else
    for v=1:3
      acc=[trjtemp.(['ka' xyz{v}])];
      accData(i,0+v)=nanmean(acc);
      accData(i,3+v)=nanstd(acc);
      accData(i,6+v)=kurtosis(acc);
    end
  end
end

%% plotting


figure(5487897)

subplot(3,2,1)
hold off
plot(velFil,velData(:,1:3),'x')
ylabel('mean velocity')
legend({'<u_x>','<u_y>','<u_z>'})

subplot(3,2,2)
hold off
plot(accFil,accData(:,1:3),'x')
ylabel('mean acceleration')
legend({'<a_x>','<a_y>','<a_z>'})

subplot(3,2,3)
hold off
plot(velFil,velData(:,4:6),'x')
ylabel('std(velocity)')
legend({'u_x''','u_y''','u_z'''})

subplot(3,2,4)
hold off
plot(accFil,accData(:,4:6),'x')
ylabel('std(acceleration)')
legend({'std(a_x)','std(a_y)','std(a_z)'})

subplot(3,2,5)
hold off
plot(velFil,velData(:,7:9),'x')
xlabel('Filter Width [frames]')
ylabel('Flatness(velocity)')
legend({'F(u_x)','F(u_y)','F(u_z)'})


subplot(3,2,6)
hold off
plot(accFil,accData(:,7:9),'x')
xlabel('Filter Width [frames]')
ylabel('Flatness(acceleration)')
legend({'F(a_x)','F(a_y)','F(a_z)'})


for k=1:6
  subplot(3,2,k)
  grid on
  set(gca,'XMinorTick','on','YMinorTick','on')
end



result.acc=accData;
result.filterAcc=accFil;
result.vel=velData;
result.filterVel=velFil;

end

