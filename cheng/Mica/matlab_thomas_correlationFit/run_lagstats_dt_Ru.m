%% Compute autocorrelation and structure functions
dtmax=100; quantity='vel';
[Rx,Sx,Npts,Ntrack]=lagstats_tracks_dt(ntraj,'x',2,dtmax,quantity,'t');
[Ry,Sy,~,~]=lagstats_tracks_dt(ntraj,'y',2,dtmax,quantity,'t');
[Rz,Sz,~,~]=lagstats_tracks_dt(ntraj,'z',2,dtmax,quantity,'t');



%% Variables
kmin=1; kmax=1000; dk=1;
% kmin=1; kmax=50; dk=5;
imin=1; imax=20;
error_max=0.002; subset_size=7;
% mRx=Rx; mRy=Ry; mRz=Rz;

%% Ru
mRx=dt1_fit(Rx,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mRy=dt1_fit(Ry,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mRz=dt1_fit(Rz,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
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

%% sv
svx=sqrt(mRx(1))*1e3; svy=sqrt(mRy(1))*1e3; svz=sqrt(mRz(1))*1e3;

%% From S to Ru
mRx_S=1-mSx/(2*mRx(kmin));
mRy_S=1-mSy/(2*mRy(kmin));
mRz_S=1-mSz/(2*mRz(kmin));
%%
t=(0:length(mRx_S)-1)/fps;
figure;hold on;grid on
%ylim([0,1])
plot(t,mRx_S);plot(t,mRy_S);plot(t,mRz_S)



%% n factors
mn21x=dt1_fit(n21x,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mn22x=dt1_fit(n22x,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mn21y=dt1_fit(n21y,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mn22y=dt1_fit(n22y,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mn21z=dt1_fit(n21z,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mn22z=dt1_fit(n22z,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
%%
color_map=lines(3);
figure;hold on;grid on
plot(sqrt(mn21x),'Color',color_map(1,:));plot(sqrt(mn22x),'Color',color_map(1,:));plot((sqrt(mn21x)+sqrt(mn22x))/2,'Color',color_map(1,:))
plot(sqrt(mn21y),'Color',color_map(2,:));plot(sqrt(mn22y),'Color',color_map(2,:));plot((sqrt(mn21y)+sqrt(mn22y))/2,'Color',color_map(2,:))
plot(sqrt(mn21z),'Color',color_map(3,:));plot(sqrt(mn22z),'Color',color_map(3,:));plot((sqrt(mn21z)+sqrt(mn22z))/2,'Color',color_map(3,:))

%% nRu
nRx=mRx./((mn21x.*mn22x).^0.5*1e-6);
nRy=mRy./((mn21y.*mn22y).^0.5*1e-6);
nRz=mRz./((mn21z.*mn22z).^0.5*1e-6);
%%
t=(0:length(nRx)-1)/fps;
figure;hold on;grid on
ylim([0,1])
plot(t,nRx);plot(t,nRy);plot(t,nRz)

%% From nS to nRu
nRx_S=1-nSx/(2*mn21x(kmin)*1e-6);
nRy_S=1-nSy/(2*mn21y(kmin)*1e-6);
nRz_S=1-nSz/(2*mn21z(kmin)*1e-6);
%%
t=(0:length(nRx_S)-1)/fps;
figure;hold on;grid on
ylim([0,1])
plot(t,nRx_S);plot(t,nRy_S);plot(t,nRz_S)

%% Fit
%two layers
%eqn='(t1*exp(-x/t1)-t2*exp(-x/t2))/(t1-t2)';
%infinite layers
eqn='1/(2*erfc(t2/t1)*exp(abs(x)/t1))*(1+erf(abs(x)/(2*t2)-t2/t1)+exp(2*abs(x)/t1)*erfc(abs(x)/(2*t2)+t2/t1))';
ft=fittype(eqn);
startpoints=[0.15 0.01];lowerpoints=[0 0];upperpoints=[1 1];
xmin=kmin; xmax=100;
color_map=lines(3);
figure;hold on;grid on;box on
set(gca,'FontSize',30,'defaultLineLineWidth',3,'defaultLineMarkerSize',8)
ylim([0,1])
xlabel('$t$(s)');ylabel('$\mathcal{R}_{uu}$')

y=nRx; x=(0:length(y)-1)'./fps;
x=x(xmin:xmax); y=y(xmin:xmax);
xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
coeff=coeffvalues(xyfit);
y=nRx; x=(0:length(y)-1)'./fps;
yfit=feval(xyfit,x);
plot(x,y,'-',x,yfit,'--','Color',color_map(1,:))
%plot(y,'-','Color',color_map(1,:));plot(yfit,'--','Color',color_map(1,:))
t1x=coeff(1)*1e3; t2x=coeff(2)*1e3;

y=nRy; x=(0:length(y)-1)'./fps;
x=x(xmin:xmax); y=y(xmin:xmax);
xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
coeff=coeffvalues(xyfit);
y=nRy; x=(0:length(y)-1)'./fps;
yfit=feval(xyfit,x);
plot(x,y,'-',x,yfit,'--','Color',color_map(2,:))
%plot(y,'-','Color',color_map(2,:));plot(yfit,'--','Color',color_map(2,:))
t1y=coeff(1)*1e3; t2y=coeff(2)*1e3;

y=nRz; x=(0:length(y)-1)'./fps;
x=x(xmin:xmax); y=y(xmin:xmax);
xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
coeff=coeffvalues(xyfit);
y=nRz; x=(0:length(y)-1)'./fps;
yfit=feval(xyfit,x);
plot(x,y,'-',x,yfit,'--','Color',color_map(3,:))
%plot(y,'-','Color',color_map(3,:));plot(yfit,'--','Color',color_map(3,:))
t1z=coeff(1)*1e3; t2z=coeff(2)*1e3;
