function [pivf,meanvel,stdvel,mtauwxy,Rv,Sf,Es] = lemanalysis(name)
%UNTITLED2 Summary of this function goes here
%   [pivf,meanvel,stdvel,mtauwxy,Rv,Sf] = lemanalysis( )
load (name)
pivf=medianf(piv);
clear piv;
nu=1e-6;

%% Extract vx and vy fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SS=size(pivf(1).vx);
N=numel(pivf); 
vy=zeros([SS N]);
vx=vy;
for k=1:N;vx(:,:,k)=pivf(k).vx;end;
for k=1:N;vy(:,:,k)=pivf(k).vy;end;
 
%% Compute average and rms vector fileds and isotropy rate %%%%%%%%%%%%%%%%
mvx=mean(vx,3);
mvy=mean(vy,3);
svx=std(vx,[],3);
svy=std(vy,[],3);
% Gives the mean of the vector field
meanvel=pivf(1);
meanvel.vx=mvx;
meanvel.vy=mvy;
meanvel.name='';
meanvel.setname='Mean Vector Field';
figure;showvec(meanvel)
% Gives the standard deviation of the vector field
stdvel=pivf(1);
stdvel.vx=svx;
stdvel.vy=svy;
stdvel.name='';
stdvel.setname='Standar Deviation Vector Field';
figure;showvec(stdvel)
%
tauxy=svx./svy;
II=find(abs(tauxy-1)<.8);
mtauwxy=mean(tauxy(II));

%% Compute correlation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rv.Rvxx=corrf(vec2scal(pivf,'ux'),1);
Rv.Rvyy=corrf(vec2scal(pivf,'uy'),2);
Rv.Rvxy=corrf(vec2scal(pivf,'ux'),2);
Rv.Rvyx=corrf(vec2scal(pivf,'uy'),1);
%% graphic corfun
figure; hold on;
plot (Rv.Rvxx.r(2:201), Rv.Rvxx.f(2:201), 'ko')
plot (Rv.Rvyy.r(2:151), Rv.Rvyy.f(2:151), 'bo')
plot (Rv.Rvxy.r(2:151), Rv.Rvxy.f(2:151), 'ro')
plot (Rv.Rvyx.r(2:201), Rv.Rvyx.f(2:201), 'go')
% fit corfun
[prxx,srxx]=polyfit(Rv.Rvxx.r(2:201), Rv.Rvxx.f(2:201),5);
frxx=polyval(prxx,Rv.Rvxx.r(2:201));
[pryy,sryy]=polyfit(Rv.Rvyy.r(2:151), Rv.Rvyy.f(2:151),5);
fryy=polyval(pryy,Rv.Rvyy.r(2:151));
[prxy,srxy]=polyfit(Rv.Rvxy.r(2:151), Rv.Rvxy.f(2:151),5);
frxy=polyval(prxy,Rv.Rvxy.r(2:151));
[pryx,sryx]=polyfit(Rv.Rvyx.r(2:201), Rv.Rvyx.f(2:201),5);
fryx=polyval(pryx,Rv.Rvyx.r(2:201));
% graphic corfun
plot (Rv.Rvxx.r(2:201), frxx, 'k-')
plot (Rv.Rvyy.r(2:151), fryy, 'b-')
plot (Rv.Rvxy.r(2:151), frxy, 'r-')
plot (Rv.Rvyx.r(2:201), fryx, 'g-')

%% Compute scalar structure function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sf.Sxx=ssf(vec2scal(pivf,'ux'),1);
Sf.Syy=ssf(vec2scal(pivf,'uy'),2);
Sf.Sxy=ssf(vec2scal(pivf,'ux'),2);
Sf.Syx=ssf(vec2scal(pivf,'uy'),1); 
%%
figure; hold on;
loglog(Sf.Sxx.r*Sf.Sxx.scaler, Sf.Sxx.sfabs(:,1),'ko')
loglog(Sf.Syy.r*Sf.Syy.scaler, Sf.Syy.sfabs(:,1),'bo')
loglog(Sf.Sxy.r*Sf.Sxy.scaler, Sf.Sxy.sfabs(:,1),'ro')
loglog(Sf.Syx.r*Sf.Syx.scaler, Sf.Syx.sfabs(:,1),'go')
%
x=.2:.01:8;
y=exp((2/3)*log(x)-4.5);
loglog(x,y,'k')

%% Calculo de parametro turbulencia 
sr23.xx=Sf.Sxx.sfabs(:,1)./((Sf.Sxx.r.*Sf.Sxx.scaler*1e-3).^(2/3))';
figure; hold on;
plot (Sf.Sxx.r*Sf.Sxx.scaler, sr23.xx, 'k')
%
sr23.yy=Sf.Syy.sfabs(:,1)./((Sf.Syy.r.*Sf.Syy.scaler*1e-3).^(2/3))';
plot (Sf.Syy.r*Sf.Syy.scaler, sr23.yy, 'b')
%
sr23.xy=Sf.Sxy.sfabs(:,1)./((Sf.Sxy.r.*Sf.Sxy.scaler*1e-3).^(2/3))';
plot (Sf.Sxy.r*Sf.Sxy.scaler, sr23.xy, 'r')
%
sr23.yx=Sf.Syx.sfabs(:,1)./((Sf.Syx.r.*Sf.Syx.scaler*1e-3).^(2/3))';
plot (Sf.Syx.r*Sf.Syx.scaler, sr23.yx, 'g')
%
mxsr23.xx=max(sr23.xx);
mxsr23.yy=max(sr23.yy);
mxsr23.xy=max(sr23.xy);
mxsr23.yx=max(sr23.yx);
%
eps.xx=(mxsr23.xx/2.1)^1.5;
eps.yy=(mxsr23.yy/2.1)^1.5;
eps.xy=(mxsr23.xy/2.1)^1.5;
eps.yx=(mxsr23.yx/2.1)^1.5;
%
eta.xx=(nu^3/eps.xx)^.25;
eta.yy=(nu^3/eps.yy)^.25;
eta.xy=(nu^3/eps.xy)^.25;
eta.yx=(nu^3/eps.yx)^.25;
%
tau_eta.xx=(nu/eps.xx)^.5;
tau_eta.yy=(nu/eps.yy)^.5;
tau_eta.xy=(nu/eps.xy)^.5;
tau_eta.yx=(nu/eps.yx)^.5;

%%
s3.xx=Sf.Sxx.sf(:,3)./((Sf.Sxx.r.*Sf.Sxx.scaler*1e-3))';
figure;plot (Sf.Sxx.r*Sf.Sxx.scaler, -s3.xx,'g')

%%
figure;plot(Sf.Sxy.sfabs(:,2)./Sf.Sxx.sfabs(:,2));

%% Power spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Es.Esx=specf(vec2scal(pivf,'ux'));
Es.Esy=specf(vec2scal(pivf,'uy'));
%%
figure;hold on;
loglog(Es.Esx.kx, Es.Esx.ex,'b');
loglog(Es.Esy.ky, Es.Esy.ey,'r');
%
x=.11:.01:4;
y=exp((-5/3)*log(x)-8.6);
loglog(x,y,'k')

end

