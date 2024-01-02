%% Compute autocorrelation and structure functions
dtmax=100; quantity='vel';
[Rx,Sx,Npts,Ntrack]=lagstats_tracks_dt(ntraj,'x',2,dtmax,quantity,'t');
[Ry,Sy,~,~]=lagstats_tracks_dt(ntraj,'y',2,dtmax,quantity,'t');
[Rz,Sz,~,~]=lagstats_tracks_dt(ntraj,'z',2,dtmax,quantity,'t');



%% Variables
kmin=3; kmax=1000; dk=1;
% kmin=1; kmax=50; dk=5;
imin=1; imax=20;
error_max=0.002; subset_size=7;
% mSx=Sx; mSy=Sy; mSz=Sz;

%% S2
mSx=dt1_fit(Sx,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mSy=dt1_fit(Sy,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mSz=dt1_fit(Sz,kmin,kmax,dk,imin,imax,error_max,subset_size,fps);
mSx=mSx*1e-6; mSy=mSy*1e-6; mSz=mSz*1e-6; %from mm to m
%%
t=(0:length(mSx)-1)/fps;
figure;hold on;grid on
set(gca,'XScale','log','YScale','log')
%xlim([1e-4,1e0]);ylim([1e-2,1e4])
plot(t(kmin:end),mSx(kmin:end));plot(t(kmin:end),mSy(kmin:end));plot(t(kmin:end),mSz(kmin:end))

%% C0
epsilon=22*1e-3; tauK=(1e-6/epsilon)^0.5;
figure;hold on;grid on
set(gca,'XScale','log','YScale','log')
%C0
%plot(t,mSx'./(epsilon*t));plot(t,mSy'./(epsilon*t));plot(t,mSz'./(epsilon*t))
plot(t/tauK,mSx'./(epsilon*t));plot(t/tauK,mSy'./(epsilon*t));plot(t/tauK,mSz'./(epsilon*t))
%<a^2>
%plot(t,mSx'./t.^2);plot(t,mSy'./t.^2);plot(t,mSz'./t.^2)

%% C0*
epsilon=22*1e-3; tau0=10.24*1e-3;
figure;hold on;grid on
set(gca,'XScale','log','YScale','log')
%plot(t,movmean(gradient(mSx),10)/epsilon*fps);plot(t,movmean(gradient(mSy),10)/epsilon*fps);plot(t,movmean(gradient(mSz),10)/epsilon*fps)
plot(t/tau0,movmean(gradient(mSx),20)/epsilon*fps);plot(t/tau0,movmean(gradient(mSy),20)/epsilon*fps);plot(t/tau0,movmean(gradient(mSz),20)/epsilon*fps)



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

%% nS
nSx=2*mSx./(mn21x+mn22x)*mn21x(kmin);
nSy=2*mSy./(mn21y+mn22y)*mn21y(kmin);
nSz=2*mSz./(mn21z+mn22z)*mn21z(kmin);
%%
t=(0:length(nSx)-1)/fps;
figure;hold on;grid on
set(gca,'XScale','log','YScale','log')
plot(t,nSx);plot(t,nSy);plot(t,nSz)

%% C0
epsilon=22*1e-3; tauK=(1e-6/(epsilon))^0.5;
figure;hold on;grid on
set(gca,'XScale','log','YScale','log')
%C0
%plot(t,nSx'./(epsilon*t));plot(t,nSy'./(epsilon*t));plot(t,nSz'./(epsilon*t))
plot(t/tauK,nSx'./(epsilon*t));plot(t/tauK,nSy'./(epsilon*t));plot(t/tauK,nSz'./(epsilon*t))
%<a^2>
%plot(t,nSx'./t.^2);plot(t,nSy'./t.^2);plot(t,nSz'./t.^2)

%% C0*
epsilon=310*1e-3; tau0=3.52*1e-3;
figure;hold on;grid on
set(gca,'XScale','log','YScale','log')
%plot(t,movmean(gradient(mSx),10)/epsilon*fps);plot(t,movmean(gradient(mSy),10)/epsilon*fps);plot(t,movmean(gradient(mSz),10)/epsilon*fps)
plot(t/tau0,movmean(gradient(nSx),30)/epsilon*fps);plot(t/tau0,movmean(gradient(nSy),30)/epsilon*fps);plot(t/tau0,movmean(gradient(nSz),30)/epsilon*fps)
