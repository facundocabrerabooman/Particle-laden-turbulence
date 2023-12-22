clear all;clc;close all

%% outlines:
% 1). first process the tracers's data, which will give the mean velocity 
% field <u> and the rms of the velocity field <u^2>^0.5 as well as 
% turbulence quantities like the dissipation rate (epsilon) and 
% fluctuation of the fluid velocity (sigmau).
% 2). then process the big particle's data 
% 3). at the end we analysis both to to have the dimensionless parameters 
% like Re_lambda) and to see their interactions (slip velocity)

%% settings for plots
addpath(genpath('C:\Users\ChengWang\Desktop\SDT_EXP'));
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultLegendLocation','best');
set(groot, 'defaultAxesTitleFontWeight','bold')
set(groot,'defaultLegendAutoUpdate','off')

set(groot, 'defaultLineLineWidth',2)
set(groot,'defaultLineMarkerSize',8)
set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultLegendFontSize',12)
set(groot, 'defaultTextFontSize',15)

%% Tracer: output path
fpath = ['.' filesep 'tracers' filesep];

fout_lagr = [fpath 'Figures_Lagr'];
fout_euler = [fpath 'Figures_Euler'];
fout_euler_subsmean = [fpath 'Figures_Euler_subsMean'];

mkdir(fout_lagr)
mkdir(fout_euler)
mkdir(fout_euler_subsmean)

%% Tracer: load Lavision output
load([fpath filesep 'STB']);

%% Tracer: set frame rate
Fs = 2996;

%% Tracer: remove the frames at the beginning
frame0 = 0;
part(1:frame0-1)=[];

%% Tracer: build track structure 
track = part2track(part);

%% Tracer: finding long tracks
L = arrayfun(@(X)(numel(X.X)),track);
long_thres = 10;
Ilong = find(L>=long_thres);

%% Tracer: find proper filter width
[s(1), m(1), w]=findFilterWidth_PTV(track(Ilong),'X');
[s(2), m(2), w]=findFilterWidth_PTV(track(Ilong),'Y');
[s(3), m(3), w]=findFilterWidth_PTV(track(Ilong),'Z');

%% Tracer: find optimal filter length
figure;
yyaxis left
loglog(w,s(1).vx,'d-',Color=color3(1,:));hold on;
loglog(w,s(2).vx,'d-',Color=color3(2,:));
loglog(w,s(3).vx,'d-',Color=color3(3,:));
hold off

yyaxis right
loglog(w,s(1).ax,'o-',Color=color3(1,:));hold on;
loglog(w,s(2).ax,'o-',Color=color3(2,:));
loglog(w,s(3).ax,'o-',Color=color3(3,:));

legend('$u_x$','$u_y$','$u_z$','$a_x$','$a_y$','$a_z$');
% title('$\sigma(w)$')
xlabel('$w$')
yyaxis left
ylabel('$\sigma_{v}$')
yyaxis right
ylabel('$\sigma_{a}$')
grid on
axis padded

%% Tracer: set values for gaussian filter
wopt = 4;
lopt = 20;

%% Tracer: plot std
plot([wopt wopt],ylim,'--',Color=color1)
savefig_custom(fout_lagr,'Std',8,7)

%% Tracer: estimate smoothed tracks with optimal filter width
tracklong_tracer=calcVelLEM(track(Ilong),wopt,lopt,Fs);
Ine=arrayfun(@(X)(~isempty(X.Vx)),tracklong_tracer)==1;
tracklong_tracer = tracklong_tracer(Ine);
save([fpath 'tracks.mat'],'long_thres','tracklong_tracer','wopt','lopt','Fs')

save([fpath 'LagrangianStats.mat'],'s','w','wopt','lopt')

clear long_thres track Ilong s m w wopt lopt Ine part

%% Tracer: velocity and acceleration pdfs
Nbins = 256; % bins for PDF
n = 10; % spacing for PDF

disp('calculating PDFs')
pdfV(1) = mkpdf5(tracklong_tracer,'Vx',Nbins,n);
pdfV(2) = mkpdf5(tracklong_tracer,'Vy',Nbins,n);
pdfV(3) = mkpdf5(tracklong_tracer,'Vz',Nbins,n);
pdfA(1) = mkpdf5(tracklong_tracer,'Ax',Nbins,2*n);
pdfA(2) = mkpdf5(tracklong_tracer,'Ay',Nbins,2*n);
pdfA(3) = mkpdf5(tracklong_tracer,'Az',Nbins,2*n);
LagragianStats.pdfV = pdfV;
LagragianStats.pdfA = pdfA;

%% Tracer: plot PDFs
figure;
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',Color=color3(1,:));hold on;
semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',Color=color3(2,:));
semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',Color=color3(3,:));
semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'o-',Color=color3(1,:));
semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'o-',Color=color3(2,:));
semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'o-',Color=color3(3,:));

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1);

legend('$u_x$','$u_y$','$u_z$','$a_x$','$a_y$','$a_z$');
% title('$PDF$')
ylabel('$PDF$')
xlabel('$u, a$')
grid on
axis padded

% add subfigure
%
% axes('Position',[0.22 0.62 0.22 0.22]);
% semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',MarkerSize=2,Color=color3(1,:));hold on;
% semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',MarkerSize=2,Color=color3(2,:));
% semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',MarkerSize=2,Color=color3(3,:));
% 
% semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'^-',MarkerSize=2,Color=color3(1,:));
% semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'^-',MarkerSize=2,Color=color3(2,:));
% semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'^-',MarkerSize=2,Color=color3(3,:));
% 
% xpdf=linspace(-5,5,1024);
% plot(xpdf,normpdf(xpdf,0,1),'k');
% grid on
% set(gca,FontSize=12)
% xlim([-5 5])

clear xpdf

savefig_custom(fout_lagr,'PDFs',8,7)

%% Tracer: compute MSD 
disp('calculating MSDs')
MSD(1) = structFunc_struct(tracklong_tracer,'Xf',2);
MSD(2) = structFunc_struct(tracklong_tracer,'Yf',2);
MSD(3) = structFunc_struct(tracklong_tracer,'Zf',2);
LagragianStats.MSD = MSD;

%% Tracer: plot MSD
figure;
loglog(MSD(1).tau/Fs,MSD(1).mean,'d-',Color=color3(1,:));hold on
loglog(MSD(2).tau/Fs,MSD(2).mean,'d-',Color=color3(2,:));
loglog(MSD(3).tau/Fs,MSD(3).mean,'d-',Color=color3(3,:));

xMSD = linspace(1,max(MSD(1).tau),1000)/Fs;
loglog(xMSD,5e4*xMSD.^2,'--',Color=color1)

legend('x','y','z')
% title('$MSD$')
ylabel('$MSD(mm^2)$')
xlabel('$\tau(s)$')
text(1e-2,2e1,'$\tau^2$','FontSize',18,FontWeight='bold')
grid on
axis padded

clear xMSD

savefig_custom(fout_lagr,'MSDs',8,7)

%% Tracer: Lagrangian 2nd SF
disp('calculating S2L')
[S2L(1), ~, ~]= structFunc_struct(tracklong_tracer,'Vx',2);
[S2L(2), ~, ~]= structFunc_struct(tracklong_tracer,'Vy',2);
[S2L(3), ~, ~]= structFunc_struct(tracklong_tracer,'Vz',2);
LagragianStats.S2L = S2L;

%% Tracer: Lagrangian plot: S2L
figure;
loglog(S2L(1).tau/Fs,S2L(1).mean,'d-',Color=color3(1,:));hold on
loglog(S2L(2).tau/Fs,S2L(2).mean,'d-',Color=color3(2,:));
loglog(S2L(3).tau/Fs,S2L(3).mean,'d-',Color=color3(3,:));

xS2L = linspace(1,4,100)/Fs;
loglog(xS2L,2e8*xS2L.^2,'--',Color=color1)

xS2L = linspace(4,100,100)/Fs;
loglog(xS2L,3e5*xS2L.^1,'--',Color=color1)

legend('x','y','z')
% title('$S_2^L$')
ylabel('$S_2^L$')
xlabel('$\tau/s$')
text(6e-4,2e2,'$\tau^2$',FontSize=18,FontWeight='bold')
text(7e-3,4e3,'$\tau$',FontSize=18,FontWeight='bold')
grid on
axis padded

clear xS2L

savefig_custom(fout_lagr,'S2L',8,7)

%% Tracer: correlation function (fit)

disp('calculating Correlation functions')
Ruu(1) = xcorr_struct(tracklong_tracer,'Vx',1);
Ruu(2) = xcorr_struct(tracklong_tracer,'Vy',1);
Ruu(3) = xcorr_struct(tracklong_tracer,'Vz',1);
Raa(1) = xcorr_struct(tracklong_tracer,'Ax',1);
Raa(2) = xcorr_struct(tracklong_tracer,'Ay',1);
Raa(3) = xcorr_struct(tracklong_tracer,'Az',1);

LagragianStats.Ruu = Ruu;
LagragianStats.Raa = Raa;

%% Tracer: fit correlation  -- Thomas's 

nCorrFitV = 300;
nCorrFitA = 20;

% seems 2layers works better for velocity 
% while infinite layers works for accerleartion
if2layersV = 1;
if2layersA = 0;

% sometimes the fit doesn't work, set the bounded option to 0 to let it work
ifboundedV = [1 1 1]; 
ifboundedA = [1 1 1];

Ruufit(1) = correlationFit(Ruu(1),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(1));
Ruufit(2) = correlationFit(Ruu(2),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(2));
Ruufit(3) = correlationFit(Ruu(3),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(3));
Raafit(1) = correlationFit(Raa(1),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(1));
Raafit(2) = correlationFit(Raa(2),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(2));
Raafit(3) = correlationFit(Raa(3),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(3));
clear nCorrFitV nCorrFitA
clear if2layersV  if2layersA
clear ifboundedV ifboundedA

%% Tracer: plot Correlation function (fit)
% figure;
% main plot: zoom in
f1 = figure;
tiledlayout(2,1)
nexttile
plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',Color=color3(1,:),MarkerSize=1);hold on
plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',Color=color3(2,:),MarkerSize=1);
plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',Color=color3(3,:),MarkerSize=1);
plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:));hold on
plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:))
plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:))

legend('$x$','$y$','$z$',Location='eastoutside')
title('$R_{uu}$')
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$');
xlabel('$\tau$/s')
grid on
axis tight

nexttile
plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'o',Color=color3(1,:),MarkerSize=1);hold on
plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'o',Color=color3(2,:),MarkerSize=1);
plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'o',Color=color3(3,:),MarkerSize=1);
plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:))
plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:))
plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:))

legend('$x$','$y$','$z$',Location='eastoutside')
title('$R_{aa}$')
ylabel('$\frac{\langle a(t)a(t+\tau) \rangle} {\langle a^2(t) \rangle}$');
xlabel('$\tau$/s')
grid on
axis tight

% add inset: zoom out 
% axes('Position',[0.4 0.3 0.3 0.2]);
% plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',Color=color3(1,:));hold on
% plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',Color=color3(2,:));
% plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',Color=color3(3,:));
% 
% plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'^',Color=color3(1,:));
% plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'^',Color=color3(2,:));
% plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'^',Color=color3(3,:));
% 
% plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:),);hold on
% plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:))
% plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:))
% 
% plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:))
% plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:))
% plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:))
% set(gca,FontSize=12)
% grid on
% axis tight

savefig_custom(fout_lagr,'Corr',8,7)

%% Tracer: find optimal ndt for Correlation function -- dt method
nmaxdt = 10;
nmaxtau = 10;
find_optimal_ndt(tracklong_tracer,nmaxdt,nmaxtau,'X',1)
find_optimal_ndt(tracklong_tracer,nmaxdt,nmaxtau,'Y',1)
find_optimal_ndt(tracklong_tracer,nmaxdt,nmaxtau,'Z',1)
clear nmaxtau nmaxtau

%% Tracer: correlation function -- dt method (denosied)
ndts = [6 6 6]; % start points of correlation function ndt
ndtl = [10 10 10]; % length of correlation function ndt

disp('calculating Correlation functions -- dt method')
[tau,corrv,corra] = dtCorr(tracklong_tracer,ndts,ndtl,Fs);

LagragianStats.ndts = ndts;
LagragianStats.ndtl = ndtl;
LagragianStats.tau = tau;
LagragianStats.corrv = corrv;
LagragianStats.corra = corra;
clear ndts ndtl

%% Tracer: plot Correlation function -- dt method (denosied) 
fields = fieldnames(tau);
figure
subplot(2,1,1)
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau.(f),corrv.(f)/corrv.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end
legend('$x$','$y$','$z$',Location='eastoutside')
xlabel('$\tau(s)$');
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$');
grid;

subplot(2,1,2)
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau.(f),corra.(f)/corra.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end
grid;
legend('$x$','$y$','$z$',Location='eastoutside')
xlabel('$\tau(s)$');
ylabel('$\frac{\langle a(t)a(t+\tau) \rangle} {\langle a^2(t) \rangle}$');

savefig_custom(fout_lagr,'dtCorr',8,7)

%% Tracer: save Lagrangian Stats
save([fpath 'LagrangianStats.mat'],'LagragianStats','-append')
disp('---- Lagrangian Stats done ----')


%% Tracer: Eulerian statistics of RAW fields  
tic
[eulerStats,pair]= twoPointsEulerianStats_Mica_Speedup(tracklong_tracer,[0.5 80],100,0,0);
save('EulerianStats.mat','eulerStats')   
save('EulerianPair.mat','pair','-v7.3')
toc

%% Tracer: compute the mean fields and substract it
dt = [4 6 8 10];
nbins = [20 21 22];
threshold = 10;
gridRange.x = [-40 40];
gridRange.y = [-40 40];
gridRange.z = [-40 40];

[gridsV,meanFieldsV,tracklong_tracer] = meanFields(tracklong_tracer,Fs,dt,nbins,threshold,1,1,gridRange,1);
[gridsVrms,meanFieldsVrms,tracklong_tracer] = meanFields(tracklong_tracer,Fs,dt,nbins,threshold,1,2,gridRange,1);
[gridsA,meanFieldsA,tracklong_tracer] = meanFields(tracklong_tracer,Fs,dt,nbins,threshold,2,1,gridRange,1);
[gridsArms,meanFieldsArms,tracklong_tracer] = meanFields(tracklong_tracer,Fs,dt,nbins,threshold,2,2,gridRange,1);
% ------>>>>>>>>
% now the tracklong contains the rms fields

%% Tracer: visualize the mean fields
slices.x = [-20 0 20];
slices.y = [0];
slices.z = [-5];
axisrange = [-40 40 -10 40 -20 20];

sliceFields(gridsV,meanFieldsV,slices,axisrange,1,1)
sliceFields(gridsVrms,meanFieldsVrms,slices,axisrange,1,2)
sliceFields(gridsA,meanFieldsA,slices,axisrange,2,1)
sliceFields(gridsArms,meanFieldsArms,slices,axisrange,2,2)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)

%% Tracer: clear useless variables 
clearvars -except Fs tracklong_tracer gridsV gridsVrms ...
    meanFieldsV meanFieldsVrms mycolormap

%% Tracer: compute the fluctuation of fluid velocity sigma_u = sqrt(<u'^2>)
sigmau1.x = sqrt(mean(vertcat(tracklong_tracer.sVx).^2,"omitnan"));
sigmau1.y = sqrt(mean(vertcat(tracklong_tracer.sVy).^2,"omitnan"));
sigmau1.z = sqrt(mean(vertcat(tracklong_tracer.sVz).^2,"omitnan"));
sigmau1.all = sqrt((sigmau1.x^2+sigmau1.y^2+sigmau1.z^2)/3); % unit: mm/s

%% Tracer: compute the fluctuation of fluid velocity sigma_u = sqrt(<u^2>-<u>^2)
[sigmauFields,sigmau2] = rmsFlucVelo(meanFieldsV,meanFieldsVrms);
sigmau2.all = sqrt((sigmau2.x^2+sigmau2.y^2+sigmau2.z^2)/3); % unit: mm/s

%% Tracer: get the rms fields
% track_subsMean: contains only the rms fields
track_subsMean = trackSubsMean(tracklong_tracer);
save('tracks_subsMean.mat','track_subsMean')

%% Tracer: Eulerian statistics of RMS fields
tic
[eulerStats_subsMean,pair_subsMean]= twoPointsEulerianStats_Mica_Speedup(track_subsMean,[0.5 80],100,1,0);
save('EulerianStats_subsMean.mat','eulerStats_subsMean')   
save('EulerianPair_subsMean.mat','pair_subsMean','-v7.3')
toc

%% Tracer: Eulerian plots of RAW or RMS field
ifsubsMean = 0;
if ifsubsMean == 0
    load('EulerianStats.mat')
else
    load('EulerianStats_subsMean.mat')
end

%% Tracer: Eulerian plots: check stationary

figure
t=tiledlayout(4,1,'TileSpacing','tight');
nexttile;
plot(eulerStats.Vmoy,'d-',Color=color1(1,:),MarkerSize=3);hold on
plot(eulerStats.VmoyX,'-',Color=color3(1,:),MarkerSize=3);
plot(eulerStats.VmoyY,'-',Color=color3(2,:),MarkerSize=3);
plot(eulerStats.VmoyZ,'-',Color=color3(3,:),MarkerSize=3);
set(gca,FontSize=15)
legend('$\langle \sqrt{x^2+y^2+z^2} \rangle$','$\langle x \rangle$','$\langle y \rangle$','$\langle z \rangle$',FontSize=10);
% title('$|V|, \sigma_V, |A|,\sigma_A$')
ylabel('$|u|$')
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Vstd,'d-',Color=color1(1,:),MarkerSize=3);hold on
plot(eulerStats.VstdX,'-',Color=color3(1,:),MarkerSize=3);
plot(eulerStats.VstdY,'-',Color=color3(2,:),MarkerSize=3);
plot(eulerStats.VstdZ,'-',Color=color3(3,:),MarkerSize=3);
set(gca,FontSize=15)
ylabel('$\sigma_u$')
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Amoy,'d-',Color=color1(1,:),MarkerSize=3);hold on
plot(eulerStats.AmoyX,'-',Color=color3(1,:),MarkerSize=3);
plot(eulerStats.AmoyY,'-',Color=color3(2,:),MarkerSize=3);
plot(eulerStats.AmoyZ,'-',Color=color3(3,:),MarkerSize=3);
set(gca,FontSize=15)
ylabel('$|a|$')
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Astd,'d-',Color=color1(1,:),MarkerSize=3);hold on
plot(eulerStats.AstdX,'-',Color=color3(1,:),MarkerSize=3);
plot(eulerStats.AstdY,'-',Color=color3(2,:),MarkerSize=3);
plot(eulerStats.AstdZ,'-',Color=color3(3,:),MarkerSize=3);
set(gca,FontSize=15)
ylabel('$\sigma_a$')
xlabel('$t/s$')
axis tight
xticks(0:Fs/10:size(eulerStats.Astd,2))
xticklabels(num2cell([0:Fs/10:size(eulerStats.Astd,2)]/Fs))
grid on

linkaxes(t.Children,'x')

savefig_custom(fout_euler,'VAt',8,7)

%% Tracer: Eulerian plots: 2nd order structure functio (S2x, S2y, S2z)

figure;
loglog(eulerStats.r,eulerStats.S2x,'d-',Color=color3(1,:));hold on
loglog(eulerStats.r,eulerStats.S2y,'d-',Color=color3(2,:));
loglog(eulerStats.r,eulerStats.S2z,'d-',Color=color3(3,:));

loglog(eulerStats.r,7e3*eulerStats.r.^(2/3),'--',Color=color1)

legend('x','y','z')
% title('$S_2^E$')
ylabel('$S_2^E$')
xlabel('$r/mm$')
text(3,3e4,'$r^{2/3}$',FontSize=18,FontWeight='bold')
grid on
axis padded

savefig_custom(fout_euler,'S2E',8,7)

%% Tracer: Eulerian plots: Ruu 
figure;
plot(eulerStats.Ruur,eulerStats.Ruu,'-',Color=color3(1,:));hold on
set(gca,FontSize=15)
% legend('$Ruu(r)$')
% title('$Ruu(r)$')
ylabel('$Ruu(r)$')
xlabel('$Ruu_r$')
grid on
axis padded

savefig_custom(fout_euler,'Ruu',8,7)

%% Tracer: Eulerian plots: PSD
figure
loglog(eulerStats.PSDk,eulerStats.PSD,'-',Color=color3(1,:));hold on

% title('$PSD$')
ylabel('$PSD$')
xlabel('$k$')
grid on
axis padded

savefig_custom(fout_euler,'PSD',8,7)

%% Tracer: Eulerian plots: Sau
figure
semilogx(eulerStats.r,eulerStats.Sau,'d-',Color=color3(1,:));hold on
semilogx(eulerStats.r,eulerStats.Saulong,'d-',Color=color3(3,:));

legend('$S_{au}$','$S_{au}^{\parallel}$')
% title('$S_{au}, S_{au}^{\parallel}$')
ylabel('$S_{au}, S_{au}^{\parallel}$')
xlabel('$r/mm$')
grid on
axis padded

savefig_custom(fout_euler,'Sau',8,7)

%% Tracer: Eulerian plots: Splong
figure
loglog(eulerStats.r,eulerStats.Splong{1,1},'d-',Color=color5(1,:));hold on
loglog(eulerStats.r,eulerStats.Splong{1,2},'d-',Color=color5(2,:));
loglog(eulerStats.r,eulerStats.Splong{1,3},'d-',Color=color5(3,:));
loglog(eulerStats.r,eulerStats.Splong{1,4},'d-',Color=color5(4,:));
loglog(eulerStats.r,eulerStats.Splong{1,5},'d-',Color=color5(5,:));

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1)

legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$')
% title('$S_n^{\parallel}$')
ylabel('$S_n^{\parallel}$')
xlabel('$r/mm$')
text(8,5e12,'$r^{5/3}$')
text(8,2e10,'$r^{4/3}$')
text(8,5e7,'$r^{3/3}$')
text(8,3e5,'$r^{2/3}$')
text(8,5e2,'$r^{1/3}$')
grid on
axis padded

savefig_custom(fout_euler,'Splong',8,7)

%% Tracer: Eulerian plots: SplongAbs
figure;
loglog(eulerStats.r,eulerStats.SplongAbs{1,1},'d-',Color=color5(1,:));hold on
loglog(eulerStats.r,eulerStats.SplongAbs{1,2},'d-',Color=color5(2,:));
loglog(eulerStats.r,eulerStats.SplongAbs{1,3},'d-',Color=color5(3,:));
loglog(eulerStats.r,eulerStats.SplongAbs{1,4},'d-',Color=color5(4,:));
loglog(eulerStats.r,eulerStats.SplongAbs{1,5},'d-',Color=color5(5,:));

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1)

legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$')
% title('$S_n^{\parallel}$')
ylabel('$S_n^{\parallel}$')
xlabel('$r/mm$')
text(8,5e12,'$r^{5/3}$')
text(8,2e10,'$r^{4/3}$')
text(8,5e7,'$r^{3/3}$')
text(8,3e5,'$r^{2/3}$')
text(8,5e2,'$r^{1/3}$')
grid on
axis padded

savefig_custom(fout_euler,'SplongAbs',8,7)

%% Tracer: Eulerian plots: dissipation rate, Epsilon
Ckolomogrov = 2.1;

figure;
loglog(eulerStats.r,(eulerStats.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats.r','d-',Color=color3(1,:));

hold on
%figure
loglog(eulerStats.r,abs(eulerStats.Splong{1,3})./(4/5*eulerStats.r)','d-',Color=color3(2,:));
loglog(eulerStats.r,abs(eulerStats.Sau)./2,'d-',Color=color3(3,:));

legend('$(S_2^{\parallel}/C_k)^{3/2}\cdot r^{-1}$','$|S_3^{\parallel}|\cdot (4/5r)^{-1}$','$|S_{au}|/2$')
% title('$\epsilon$')
ylabel('$\epsilon$')
xlabel('$r/mm$')
grid on;axis padded

savefig_custom(fout_euler,'Epsilon',8,7)

%% !!!!!!! clear workspace !!!!!!! 
clearvars -except Fs tracklong_tracer gridsV gridsVrms ...
    meanFieldsV meanFieldsVrms 

%% Particle: output path
fpath = ['.' filesep 'particle' filesep];
fout_lagr = [fpath 'Figures_Lagr'];
mkdir(fout_lagr)

%% Particle: load Lavision output
load([fpath filesep 'STB']);

%% Particle: remove fields (Lavision output)
fields = {'Vx','Vy','Vz','V','Ax','Ay','Az','A'};
part = rmfield(part,fields);
clear fields

%% Particle: remove the frames at the beginning
frame0 = 0;
part(1:frame0-1)=[];

%% Particle: build track structure 
track = part2track(part);

%% Particle: finding long tracks
L = arrayfun(@(X)(numel(X.X)),track);
long_thres = 10;
Ilong = find(L>=long_thres);

%% Particle: find proper filter width
[s(1), m(1), w]=findFilterWidth_PTV(track(Ilong),'X');
[s(2), m(2), w]=findFilterWidth_PTV(track(Ilong),'Y');
[s(3), m(3), w]=findFilterWidth_PTV(track(Ilong),'Z');

%% Particle: find optimal filter length
figure;
yyaxis left
loglog(w,s(1).vx,'d-',Color=color3(1,:));hold on;
loglog(w,s(2).vx,'d-',Color=color3(2,:));
loglog(w,s(3).vx,'d-',Color=color3(3,:));
hold off

yyaxis right
loglog(w,s(1).ax,'o-',Color=color3(1,:));hold on;
loglog(w,s(2).ax,'o-',Color=color3(2,:));
loglog(w,s(3).ax,'o-',Color=color3(3,:));

legend('$u_x$','$u_y$','$u_z$','$a_x$','$a_y$','$a_z$');
% title('$\sigma(w)$')
xlabel('$w$')
yyaxis left
ylabel('$\sigma_{v}$')
yyaxis right
ylabel('$\sigma_{a}$')
grid on
axis padded

%% Particle: set values for gaussian filter
wopt = 4;
lopt = 20;

%% Particle: plot std
plot([wopt wopt],ylim,'--',Color=color1)
savefig_custom(fout_lagr,'Std',8,7)

%% Particle: Estimate smoothed tracks with optimal filter width
tracklong_particle=calcVelLEM(track(Ilong),wopt,lopt,Fs);
Ine=arrayfun(@(X)(~isempty(X.Vx)),tracklong_particle)==1;
tracklong_particle = tracklong_particle(Ine);
save([fpath 'tracks.mat'],'long_thres','tracklong_particle','wopt','lopt','Fs')

save([fpath 'LagrangianStats.mat'],'s','w','wopt','lopt')

clear long_thres track Ilong s m w wopt lopt Ine part

%% Particle: plot trajectory, V(t), A(t) -build 
particle_part_filtered = track2part(tracklong_particle,{'Tf','Xf','Yf','Zf','Vx','Vy','Vz','Ax','Ay','Az'},1);

%% Particle: plot trajectory, V(t), A(t) -plot 
fout_traj = [fpath './Figures_traj'];
titlestring = '$00V-0p6g$';
axisrange_traj = [-13 7 -15 5 -40 40];
axisrange_vat = [0.05 0.35];

plot_traj(particle_part_filtered,titlestring,axisrange_traj,1,fout_traj)
plot_vadiffusion(particle_part_filtered,axisrange_vat,Fs,1,fout_traj)

clear fout_traj titlestring axisrange_traj axisrange_vat

%% Particle: velocity and acceleration pdfs
Nbins = 256; % bins for PDF
n = 10; % spacing for PDF

disp('calculating PDFs')
pdfV(1) = mkpdf5(tracklong_particle,'Vx',Nbins,n);
pdfV(2) = mkpdf5(tracklong_particle,'Vy',Nbins,n);
pdfV(3) = mkpdf5(tracklong_particle,'Vz',Nbins,n);
pdfA(1) = mkpdf5(tracklong_particle,'Ax',Nbins,2*n);
pdfA(2) = mkpdf5(tracklong_particle,'Ay',Nbins,2*n);
pdfA(3) = mkpdf5(tracklong_particle,'Az',Nbins,2*n);
LagragianStats.pdfV = pdfV;
LagragianStats.pdfA = pdfA;

%% Particle: plot PDFs
figure;
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',Color=color3(1,:));hold on;
semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',Color=color3(2,:));
semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',Color=color3(3,:));
semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'o-',Color=color3(1,:));
semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'o-',Color=color3(2,:));
semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'o-',Color=color3(3,:));

xpdf=linspace(-4,4,1024);
plot(xpdf,normpdf(xpdf,0,1),'k--');

legend('$u_x$','$u_y$','$u_z$','$a_x$','$a_y$','$a_z$');
% title('$PDF$')
ylabel('$PDF$')
xlabel('$u, a$')
grid on
axis padded

% add subfigure
%
% axes('Position',[0.22 0.62 0.22 0.22]);
% semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',MarkerSize=2,Color=color3(1,:));hold on;
% semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',MarkerSize=2,Color=color3(2,:));
% semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',MarkerSize=2,Color=color3(3,:));
% 
% semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'^-',MarkerSize=2,Color=color3(1,:));
% semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'^-',MarkerSize=2,Color=color3(2,:));
% semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'^-',MarkerSize=2,Color=color3(3,:));
% 
% xpdf=linspace(-5,5,1024);
% plot(xpdf,normpdf(xpdf,0,1),'k');
% grid on
% set(gca,FontSize=12)
% xlim([-5 5])

clear xpdf

savefig_custom(fout_lagr,'PDFs',8,7)

%% Particle: compute MSD 
disp('calculating MSDs')
MSD(1) = structFunc_struct(tracklong_particle,'Xf',2);
MSD(2) = structFunc_struct(tracklong_particle,'Yf',2);
MSD(3) = structFunc_struct(tracklong_particle,'Zf',2);
LagragianStats.MSD = MSD;

%% Particle: plot MSD
figure;
loglog(MSD(1).tau/Fs,MSD(1).mean,'d-',Color=color3(1,:));hold on
loglog(MSD(2).tau/Fs,MSD(2).mean,'d-',Color=color3(2,:));
loglog(MSD(3).tau/Fs,MSD(3).mean,'d-',Color=color3(3,:));

xMSD = linspace(1,max(MSD(1).tau),1000)/Fs;
loglog(xMSD,5e4*xMSD.^2,'k--')

legend('x','y','z')
% title('$MSD$')
ylabel('$MSD(mm^2)$')
xlabel('$\tau(s)$')
text(1e-2,2e1,'$\tau^2$','FontSize',18,FontWeight='bold')
grid on
axis padded

clear xMSD

savefig_custom(fout_lagr,'MSDs',8,7)

%% Particle: Lagrangian 2nd SF
disp('calculating S2L')
[S2L(1), ~, ~]= structFunc_struct(tracklong_particle,'Vx',2);
[S2L(2), ~, ~]= structFunc_struct(tracklong_particle,'Vy',2);
[S2L(3), ~, ~]= structFunc_struct(tracklong_particle,'Vz',2);
LagragianStats.S2L = S2L;

%% Particle: Lagrangian plot: S2L
figure;
loglog(S2L(1).tau/Fs,S2L(1).mean,'d-',Color=color3(1,:));hold on
loglog(S2L(2).tau/Fs,S2L(2).mean,'d-',Color=color3(2,:));
loglog(S2L(3).tau/Fs,S2L(3).mean,'d-',Color=color3(3,:));

xS2L = linspace(1,4,100)/Fs;
loglog(xS2L,2e8*xS2L.^2,'k--')

xS2L = linspace(4,100,100)/Fs;
loglog(xS2L,3e5*xS2L.^1,'k--')

legend('x','y','z')
% title('$S_2^L$')
ylabel('$S_2^L$')
xlabel('$\tau/s$')
text(6e-4,2e2,'$\tau^2$',FontSize=18,FontWeight='bold')
text(7e-3,4e3,'$\tau$',FontSize=18,FontWeight='bold')
grid on
axis padded

clear xS2L

savefig_custom(fout_lagr,'S2L',8,7)

%% Particle: correlation function (fit)
disp('calculating Correlation functions')
Ruu(1) = xcorr_struct(tracklong_particle,'Vx',1);
Ruu(2) = xcorr_struct(tracklong_particle,'Vy',1);
Ruu(3) = xcorr_struct(tracklong_particle,'Vz',1);
Raa(1) = xcorr_struct(tracklong_particle,'Ax',1);
Raa(2) = xcorr_struct(tracklong_particle,'Ay',1);
Raa(3) = xcorr_struct(tracklong_particle,'Az',1);

LagragianStats.Ruu = Ruu;
LagragianStats.Raa = Raa;

%% Particle: fit correlation  -- Thomas's 
nCorrFitV = [300 300 300];
nCorrFitA = [21 19 21];

% seems 2layers works better for velocity 
% while infinite layers works for accerleartion
if2layersV = 1;
if2layersA = 0;

% sometimes the fit doesn't work, set the bounded option to 0 to let it work
ifboundedV = [1 1 1]; 
ifboundedA = [0 0 0];

Ruufit(1) = correlationFit(Ruu(1),Fs,1,nCorrFitV(1),'V',if2layersV,ifboundedV(1));
Ruufit(2) = correlationFit(Ruu(2),Fs,1,nCorrFitV(2),'V',if2layersV,ifboundedV(2));
Ruufit(3) = correlationFit(Ruu(3),Fs,1,nCorrFitV(3),'V',if2layersV,ifboundedV(3));
Raafit(1) = correlationFit(Raa(1),Fs,1,nCorrFitA(1),'A',if2layersA,ifboundedA(1));
Raafit(2) = correlationFit(Raa(2),Fs,1,nCorrFitA(2),'A',if2layersA,ifboundedA(2));
Raafit(3) = correlationFit(Raa(3),Fs,1,nCorrFitA(3),'A',if2layersA,ifboundedA(3));
clear nCorrFitV nCorrFitA
clear if2layersV  if2layersA
clear ifboundedV ifboundedA

%% Particle: plot Correlation function (fit)
% figure;
% main plot: zoom in
f1 = figure;
tiledlayout(2,1)
nexttile
plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',Color=color3(1,:),MarkerSize=1);hold on
plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',Color=color3(2,:),MarkerSize=1);
plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',Color=color3(3,:),MarkerSize=1);
plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:));hold on
plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:))
plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:))

legend('$x$','$y$','$z$',Location='eastoutside')
title('$R_{uu}$')
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$');
xlabel('$\tau$/s')
grid on
axis tight

nexttile
plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'o',Color=color3(1,:),MarkerSize=1);hold on
plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'o',Color=color3(2,:),MarkerSize=1);
plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'o',Color=color3(3,:),MarkerSize=1);
plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:))
plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:))
plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:))

legend('$x$','$y$','$z$',Location='eastoutside')
title('$R_{aa}$')
ylabel('$\frac{\langle a(t)a(t+\tau) \rangle} {\langle a^2(t) \rangle}$');
xlabel('$\tau$/s')
grid on
axis tight

% add inset: zoom out 
% axes('Position',[0.4 0.3 0.3 0.2]);
% plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',Color=color3(1,:));hold on
% plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',Color=color3(2,:));
% plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',Color=color3(3,:));
% 
% plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'^',Color=color3(1,:));
% plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'^',Color=color3(2,:));
% plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'^',Color=color3(3,:));
% 
% plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:),);hold on
% plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:))
% plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:))
% 
% plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:))
% plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:))
% plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:))
% set(gca,FontSize=12)
% grid on
% axis tight

savefig_custom(fout_lagr,'Corr',8,7)

%% Particle: find optimal ndt for Correlation function -- dt method
nmaxdt = 10;
nmaxtau = 10;
find_optimal_ndt(tracklong_particle,nmaxdt,nmaxtau,'X','save',fpath)
find_optimal_ndt(tracklong_particle,nmaxdt,nmaxtau,'Y','save',fpath)
find_optimal_ndt(tracklong_particle,nmaxdt,nmaxtau,'Z','save',fpath)
clear nmaxtau nmaxtau

%% Particle: correlation function -- dt method (denosied)
ndts = [6 6 6]; % start points of correlation function ndt
ndtl = [10 10 10]; % length of correlation function ndt

disp('calculating Correlation functions -- dt method')
[tau,corrv,corra] = dtCorr(tracklong_particle,ndts,ndtl,Fs);

LagragianStats.ndts = ndts;
LagragianStats.ndtl = ndtl;
LagragianStats.tau = tau;
LagragianStats.corrv = corrv;
LagragianStats.corra = corra;
clear ndts ndtl

%% Particle: plot Correlation function -- dt method (denosied) 
fields = fieldnames(tau);
figure
subplot(2,1,1)
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau.(f),corrv.(f)/corrv.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end
legend('$x$','$y$','$z$',Location='eastoutside')
xlabel('$\tau(s)$');
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$');
grid;

subplot(2,1,2)
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau.(f),corra.(f)/corra.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end
grid;
legend('$x$','$y$','$z$',Location='eastoutside')
xlabel('$\tau(s)$');
ylabel('$\frac{\langle a(t)a(t+\tau) \rangle} {\langle a^2(t) \rangle}$');

savefig_custom(fout_lagr,'dtCorr',8,7)

%% Particle: save Lagrangian Stats
save([fpath 'LagrangianStats.mat'],'LagragianStats','-append')
disp('---- Lagrangian Stats done ----')

%% !!!!!!! clear workspace !!!!!!! 
clearvars -except Fs tracklong_tracer gridsV gridsVrms ...
    meanFieldsV meanFieldsVrms tracklong_particle particle_part_filtered

%% BOTH: make tracking video
Nexp = 2;
trajlen = 100;
trajincrp = 10;
trajincrt = 1;

video3Dtraj(particle_part_filtered,tracklong_tracer,Nexp,trajlen,trajincrp,trajincrt,Rmin,Rmax);
save('neighbor.mat','neighbor_global','neighbor_layer')
clear Nexp trajlen trajincrp trajincrt Rmin Rmax

%% BOTH: experiment parameters
mp = 0.0048e-3;%kg
dp = 1.0094e-3;%m
rhof = 1000; %kg/m3
nv = 8.532e-7;%m2/s, @ 27 celsius
geff = 1; %ratio to 9.8m/s2
taup = 0.10; % determine from the plots of V(t) and Vs(t)
Vs = 0.4; % abs, determine from the plots of V(t) and Vs(t)
epsilon = 0.15; %determine from Eulerian plots of tracers
sigmau = sigmau2; % determine from the Eulerian fields, unit: mm/s

params = expParam(Fs,mp,dp,rhof,nv,geff,taup,Vs,epsilon,sigmau)
save('./expParams.mat','params');

clear mp dp rhof nv geff taup 

%% BOTH: slip velocity: find terminal state where a ~= 0
[maxTimeLength,startTime,endTime,ThresErrorTS] = findTerminalState(part,1);

%% BOTH: slip velocity: neighbor tracks at terminal state
Rmin = 0;
Rmax = 5;
numExp = numel(startTime);
neighborAll = [];neighborLayer=[];trackTSFront=[];trackTSBack=[];
for kexp = 1:numExp
    [neighborAllTemp, neighborLayerTemp] = neighbor(part,tracklong_particle,Rmin,Rmax,kexp);

    [trackTSFrontTemp, trackTSBackTemp]  = neighborAllTS(part,tracklong_particle,neighborAllTemp,startTime,endTime,kexp);

    nframe(kexp) = numel(neighborAllTemp);

    neighborAll   = [neighborAll;   neighborAllTemp];
    neighborLayer = [neighborLayer; neighborLayerTemp];
    trackTSFront  = [trackTSFront;  trackTSFrontTemp];
    trackTSBack   = [trackTSBack;   trackTSBackTemp];

    % data contained in trackTSFront and trackTSBack: 
    % coordinates (X,Y,Z,Xf,Yf,Zf) are in particle reference framework (Xf-Xp)
    % tracers' velocity and acceleration (vx,vy,vz,ax,ay,az) are also
    % minus particle's conterpart (tracers.Vx - part.Vx)
    % a rotaion matrix of each instant time was ALREADY applied to 
    % the relative coordinates, velocity, acceleration

    clear neighborAllTemp neighborLayerTemp trackTSFrontTemp trackTSBackTemp 
end

trackTSALL = mergeStructures(trackTSFront, trackTSBack);

save('./trackTS.mat','trackTSALL','trackTSFront','trackTSBack');

%% BOTH: slip velocity: INTERPOLATION: meanFields around particle
dt = [2 4]; % change it if error message show up
nbins = [3 4 5];
threshold = 0;
gridRange.x = [-Rmax Rmax];
gridRange.y = [-Rmax Rmax];
gridRange.z = [-Rmax Rmax];
[gridsSlipVeloInterp,slipVeloInterp,trackTSALL] = meanFields(trackTSALL,Fs,dt,nbins,threshold,1,1,gridRange,0);
% [gridsSlipAcceInterp,slipAcceInterp,trackTSALL] = meanFields(trackTSALL,Fs,dt,nbins,threshold,2,1,gridRange,0);

save('./slipVeloINTERP.mat','gridsSlipVeloInterp','slipVeloInterp');

%% BOTH: slip velocity: INTERPOLATION: visulaize the flow around particle
slices.x = [0];
slices.y = [0];
slices.z = [0];
axisrange = [-Rmax Rmax -Rmax Rmax -Rmax Rmax];

sliceFields(gridsSlipVeloInterp,slipVeloInterp,slices,axisrange,1,1,1)
% sliceFields(gridsSlipAcceInterp,slipAcceInterp,slices,axisrange,2,1,1)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)

%% BOTH: slip velocity: INTERPOLATION: slice of mean surrounding flow
nSlice = 10;
rp = 0.5; % mm
gridSpacing = 0.1; % mm
nSpacing = 2*Rmax/gridSpacing+1;
[myMeanSliceInterp,sliceDataInterp] = sliceSlipVelo(gridsSlipVeloInterp,slipVeloInterp,nSlice,rp,Rmax,nSpacing,1,1,1);

save('./slipVeloINTERP.mat','myMeanSliceInterp','sliceDataInterp','-append');

%% BOTH: slip velocity: INTERPOLATION: plot mean slice
surf(myMeanSliceInterp);
view(90,0)

%% BOTH: slip velocity: fit into BINS
dt = [2 4]; % change it if error message show up
nbins = [3 4 5];

[gridsSlipVeloBin,slipVeloBin] = slipVeloBins(trackTSALL,Fs,dt,nbins,1,1);
% [gridsSlipAcceBins,slipAcceBins] = slipVeloBins(trackTSALL,Fs,dt,nbins,2,1);

save('./slipVeloBIN.mat','gridsSlipVeloBin','slipVeloBin');

%% BOTH: slip velocity: fit into BINS: visulaize the flow around particle
slices.x = [0];
slices.y = [0];
slices.z = [0];
axisrange = [-Rmax Rmax -Rmax Rmax -Rmax Rmax];

sliceFields(gridsSlipVeloBin,slipVeloBin,slices,axisrange,1,1,1)
% sliceFields(gridsSlipAcceInterp,slipAcceInterp,slices,axisrange,2,1,1)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)

%% BOTH: slip velocity: fit into BINS: slice of mean surrounding flow
nSlice = 10;
rp = 0.5; % mm
gridSpacing = 0.1; % mm
nSpacing = 2*Rmax/gridSpacing+1;
[myMeanSliceBin,sliceDataBin] = sliceSlipVelo(gridsSlipVeloBin,slipVeloBin,nSlice,rp,Rmax,nSpacing,1,1,1);

save('./slipVeloBIN.mat','myMeanSliceBin','sliceDataBin','-append');

%% BOTH: slip velocity: fit into BINS: plot mean slice
surf(myMeanSliceBin);
view(90,0)