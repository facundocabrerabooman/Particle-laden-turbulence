%% Compute autocorrelation and structure functions
dtmax=100; quantity='acc';
[Rx,Sx,Npts,Ntrack]=lagstats_tracks_dt(ntraj,'x',2,dtmax,quantity,'t');
[Ry,Sy,~,~]=lagstats_tracks_dt(ntraj,'y',2,dtmax,quantity,'t');
[Rz,Sz,~,~]=lagstats_tracks_dt(ntraj,'z',2,dtmax,quantity,'t');



%% Variables
kmin=1; kmax=200; dk=1;
% kmin=1; kmax=50; dk=5;
imin=1; imax=11;
error_max=0.002; subset_size=3;
% mRx=Rx; mRy=Ry; mRz=Rz;

%% Ra
mRx=dt2_fit(Rx,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mRy=dt2_fit(Ry,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mRz=dt2_fit(Rz,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mRx=mRx*1e-6; mRy=mRy*1e-6; mRz=mRz*1e-6; %from mm to m
%%
t=(0:length(mRx)-1)/fps;
mR=(mRx+mRy+mRz)/3;
figure;hold on;grid on
%xlim([0,0.3]);ylim([-0.3e5,1.5e5])
%plot(mRx(kmin:end));plot(mRy(kmin:end));plot(mRz(kmin:end))
%plot(mR(kmin:end))
plot(t(kmin:end),mRx(kmin:end));plot(t(kmin:end),mRy(kmin:end));plot(t(kmin:end),mRz(kmin:end))
%plot(t(kmin:end),mR(kmin:end))

%% a0
sa2=mR(kmin);
sax2=mRx(kmin); say2=mRy(kmin); saz2=mRz(kmin);

%% t0
t0=(find(mR<0,1)-1)/fps*1e3;
t0x=(find(mRx<0,1)-1)/fps*1e3; t0y=(find(mRy<0,1)-1)/fps*1e3; t0z=(find(mRz<0,1)-1)/fps*1e3;

%% Fit
%two layers
%eqn='a*(t1*exp(-x/t2)-t2*exp(-x/t1))/(t1-t2)';
%infinite layers
eqn='a/(2*t1/(sqrt(pi)*t2)*exp(-(t2/t1)^2)-2*erfc(t2/t1))*(2*t1/(sqrt(pi)*t2)*exp(-((x/(2*t2))^2+(t2/t1)^2))-exp(-abs(x)/t1)*(1+erf(abs(x)/(2*t2)-t2/t1)) -exp(abs(x)/t1)*erfc(abs(x)/(2*t2)+t2/t1))';
ft=fittype(eqn);
startpoints=[max(mRx) 0.1 0.01]; lowerpoints=[0.9*max(mRx) 0 0]; upperpoints=[1.1*max(mRx) 0.5 0.5];
y=mR; x=(0:length(y)-1)'./fps;
xmin=15; xmax=66;
x=x(xmin:xmax); y=y(xmin:xmax);
xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
coeff=coeffvalues(xyfit);
errors=confint(xyfit);
sa2=coeff(1); sa2_error=sa2-errors(1,1);
t1=coeff(2)*1e3; t1_error=t1-errors(1,2)*1e3; %in ms
t2=coeff(3)*1e3; t2_error=t2-errors(1,3)*1e3; %in ms
y=mR; x=(0:length(y)-1)'./fps;
yfit=feval(xyfit,x);
figure;hold on;grid on;box on
set(gca,'FontSize',30,'defaultLineLineWidth',3,'defaultLineMarkerSize',8)
plot(x,y,'-',x,yfit,'--')
%plot(y,'-');plot(yfit,'--')
xlabel('$t$(s)');ylabel('$R_{aa}$(m$^2$/s$^4$)')

%old
%color_map=lines(3);
% y=mRx; x=(0:length(y)-1)'./fps;
% x=x(kmin:kmax); y=y(kmin:kmax);
% xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
% coeff=coeffvalues(xyfit);
% y=mRx; x=(0:length(y)-1)'./fps;
% yfit=feval(xyfit,x);
% plot(x,y,'-',x,yfit,'--','Color',color_map(1,:))
% sax2=coeff(1); t1x=coeff(2)*1e3; t2x=coeff(3)*1e3;
% 
% y=mRy; x=(0:length(y)-1)'./fps;
% x=x(kmin:kmax); y=y(kmin:kmax);
% xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
% coeff=coeffvalues(xyfit);
% y=mRy; x=(0:length(y)-1)'./fps;
% yfit=feval(xyfit,x);
% plot(x,y,'-',x,yfit,'--','Color',color_map(2,:))
% say2=coeff(1); t1y=coeff(2)*1e3; t2y=coeff(3)*1e3;
% 
% y=mRz; x=(0:length(y)-1)'./fps;
% x=x(kmin:kmax); y=y(kmin:kmax);
% xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
% coeff=coeffvalues(xyfit);
% y=mRz; x=(0:length(y)-1)'./fps;
% yfit=feval(xyfit,x);
% plot(x,y,'-',x,yfit,'--','Color',color_map(3,:))
% saz2=coeff(1); t1z=coeff(2)*1e3; t2z=coeff(3)*1e3;
