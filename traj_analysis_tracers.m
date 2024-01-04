clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Volumes/landau1/Tracers/ddt_moredata';
folderout = folderin;
cd(folderin)

Fs=2990; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Concatenate data

if pi==pi
    trajs_conc = [];
    
    load('trajs_TrCer_1000_11_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    1
    load('trajs_TrCer_1000_13_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    2
    load('trajs_TrCer_1000_14_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    3
    load('trajs_TrCer_1000_15_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    4
    load('trajsf1_TrCer_1000_32_ddt_tracers.mat','tracklong1')
    trajs_conc = [trajs_conc tracklong1];
    clear tracklong1
    5
    load('trajsf2_TrCer_1000_32_ddt_tracers.mat','tracklong2')
    trajs_conc = [trajs_conc tracklong2];
    clear tracklong2
    6
    load('trajsf_TrCer_1000_31_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_30_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_28_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_29_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_27_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7

    
    % load('trajs_TrCer_1000_noturb.mat')
    % trajs_conc = [trajs_conc tracklong];
    
    % load('ddt_july3.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 3
    % load('ddt_july7a.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 4
    % load('ddt_july7b.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 5
    % load('ddt_july7c.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 6
    % load('ddt_july9b.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 7

    Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
end
trajs_conc = trajs_conc(Ine);
clear tracklong traj_ddt
try
save('traj_conc_tracers_ddt','trajs_conc','-v7.3')
catch
end

%% Get Mean Velocity 
if pi==pi
for i=1:numel(trajs_conc)
    trajs_conc(i).X = trajs_conc(i).x;
    trajs_conc(i).Y = trajs_conc(i).y;
    trajs_conc(i).Z = trajs_conc(i).z;
end

[U, mBdt, bins] = track2meanDxDt3DProfile(trajs_conc,'Xf',2:2:50,[20 20 20],1,1,'x','cart');
[V, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Yf',2:2:50,[20 20 20],1,1,'y','cart');
[W, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Zf',2:2:50,[20 20 20],1,1,'z','cart');

[X,Y,Z]=meshgrid(bins{1},bins{2},bins{3});

U=U.*Fs;
V=V.*Fs;
W=W.*Fs;

%save('output_Vel_meanfields','U','V','W','mBdt','bins','X','Y','Z')


trajs_conc_minus_mean_field = find_closest_bin(trajs_conc, X, Y, Z, U, V, W);
trajs_conc = trajs_conc_minus_mean_field;

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);

trajs_conc = trajs_conc(Ine);

try
save('trajs_conc_minus_mean_field','trajs_conc_minus_mean_field','-v7.3')
catch
end
end
%% 1 time - 1 particle statistics
if pi==pi
%% Calculate & plot velocity and acceleration pdfs
pdfV(1) = mkpdf5(trajs_conc(Ine),'Vx',256,10);
1
pdfV(2) = mkpdf5(trajs_conc(Ine),'Vy',256,10);
2
pdfV(3) = mkpdf5(trajs_conc(Ine),'Vz',256,10);
3

pdfA(1) = mkpdf5(trajs_conc(Ine),'Ax',256,20);
4
pdfA(2) = mkpdf5(trajs_conc(Ine),'Ay',256,20);
5
pdfA(3) = mkpdf5(trajs_conc(Ine),'Az',256,20);
6

try
save('output_post_processing.mat','pdfV','pdfA')
catch end
%% Plot Normalized PDFs
figure;
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on;
semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);
semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d',MarkerSize=3,Color=color3(3,:),LineWidth=2);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=3);

%set(gca,FontSize=15)
legend('$V_x$','$V_y$','$V_z$','Gaussian','interpreter','latex',Location='northeast');
%title('$PDF$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
xlabel('$\frac{V-\langle V \rangle}{std(V)}$','interpreter','latex',FontWeight='bold')
grid off
axis padded
ylim([5e-7 1])

% text(5,1,['MeanAX = ' num2str(pdfA(1).mean)])
% text(5,0.6,['MeanAY = ' num2str(pdfA(2).mean)])
% text(5,0.3,['MeanAZ = ' num2str(pdfA(3).mean)])
% 
% text(5,0.1,['MeanVX = ' num2str(pdfV(1).mean)])
% text(5,0.05,['MeanVY = ' num2str(pdfV(2).mean)])
% text(5,0.03,['MeanVZ = ' num2str(pdfV(3).mean)])

folderout = 'pdfs/';
mkdir(folderout)
savefig_FC([folderout filesep 'PDF_v'],8,6,'pdf')
savefig_FC([folderout filesep 'PDF_v'],8,6,'fig')
%%%%%%%%%

figure;

semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2); hold on
semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);
semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'d',MarkerSize=3,Color=color3(3,:),LineWidth=2);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=3);

ylabel('$PDF(\frac{A-\langle A \rangle}{std(A)})$','interpreter','latex',FontWeight='bold')
legend('$A_x$','$A_y$','$A_z$','Gaussian','interpreter','latex',Location='northeast');
xlabel('$\frac{A-\langle A \rangle}{std(A)}$','interpreter','latex',FontWeight='bold')
grid off
%xlim([-5 5])
ylim([5e-7 1])


folderout = 'pdfs/';
mkdir(folderout)
savefig_FC([folderout 'PDF_a'],8,6,'pdf')
savefig_FC([folderout 'PDF_a'],8,6,'fig')

%% Table with moments of distribution
maketable(pdfA,pdfV,folderout)

%%






%%%%%%%%%%%%%%%%%% 2 times - 1 particle statistics (Lagrangian statistics)







%% Mean Square Separation
MSD(1) = structFunc_struct(trajs_conc,'Xf',2);
1
MSD(2) = structFunc_struct(trajs_conc,'Yf',2);
2
MSD(3) = structFunc_struct(trajs_conc,'Zf',2);
3

try
save('output_post_processing.mat','MSD','-append')
catch end
%%
figure;
loglog(MSD(1).tau/Fs,MSD(1).mean,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on
loglog(MSD(2).tau/Fs,MSD(2).mean,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);
loglog(MSD(3).tau/Fs,MSD(3).mean,'d',MarkerSize=3,Color=color3(3,:),LineWidth=2);

xMSD = linspace(1,350,1000)/Fs;
loglog(xMSD,0.5e5*xMSD.^2,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$MSD^x$','$MSD^y$','$MSD^z$','Interpreter','latex', 'Location','southeast')
%title('$MSD$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$MSD(m^2)$','interpreter','latex',FontWeight='bold')
xlabel('$\tau(s)$','interpreter','latex',FontWeight='bold')
text(2e-3,3,'$\tau^2$','interpreter','latex',FontWeight='bold',FontSize=20)
grid on
axis padded


folderout = 'MSS/';
mkdir(folderout)
savefig_FC([folderout 'MSS'],8,6,'pdf')
savefig_FC([folderout 'MSS'],8,6,'fig')

%% Longitudinal S2

S2L(1)= structFunc_struct(trajs_conc,'Vx',2);
S2L(2)= structFunc_struct(trajs_conc,'Vy',2);
S2L(3)= structFunc_struct(trajs_conc,'Vz',2);

try
save('output_post_processing.mat','S2L')
catch end
%%
% figure;loglog(S2Lx.tau,S2Lx.mean./S2Lx.tau/Fs/2)
figure;
loglog(S2L(1).tau/Fs,S2L(1).mean,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on
loglog(S2L(2).tau/Fs,S2L(2).mean,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);
loglog(S2L(3).tau/Fs,S2L(3).mean,'d',MarkerSize=3,Color=color3(3,:),LineWidth=2);

xS2L = linspace(1,15,100)/Fs;
loglog(xS2L,2e8*xS2L.^2,'--',Color=color1,LineWidth=2)
xS2L = linspace(16,200,100)/Fs;
loglog(xS2L,9.5e5*xS2L.^1,'--',Color=color1,LineWidth=2)
% xS2L = linspace(100,300,100);
% loglog(xS2L,8e4*xS2L.^0,'--',Color=color1,LineWidth=2)

legend('$S_2^L(x)$','$S_2^L(y)$','$S_2^L(z)$','interpreter','latex',Location='southeast')
ylabel('$S_2^L$','interpreter','latex',FontWeight='bold')
xlabel('$\tau / sec$','interpreter','latex',FontWeight='bold')
text(8e-4,8e2,'$\tau^2$','interpreter','latex','FontWeight','bold','FontSize',20)
text(1e-2,1.5e4,'$\tau$','interpreter','latex','FontWeight','bold','FontSize',20)
grid on
axis padded

folderout = 'S2L/';
mkdir(folderout)
savefig_FC([folderout 'S2L'],8,6,'pdf')
savefig_FC([folderout 'S2L'],8,6,'fig')
end
%% Velocity and Acceleration Correlations
if pi==pi
disp('corr')

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);


Ruu(1) = xcorr_struct(trajs_conc,'Vx',1);
Ruu(2) = xcorr_struct(trajs_conc,'Vy',1);
Ruu(3) = xcorr_struct(trajs_conc,'Vz',1);

Raa(1) = xcorr_struct(trajs_conc,'Ax',1);
%load('output_post_processing.mat','Ruu','Raa')
Raa(2) = xcorr_struct(trajs_conc,'Ay',1);
Raa(3) = xcorr_struct(trajs_conc,'Az',1);

try
save('output_post_processing.mat','Ruu','Raa')
catch end
%% Tracer: fit correlation  -- Thomas's 
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/cheng'));

nCorrFitV = 270;
nCorrFitA = 20;

% seems 2layers works better for velocity 
% while infinite layers works for accerleartion
if2layersV = 2;
if2layersA = 99;

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

folderout = 'corr/';
mkdir(folderout)
savefig_FC([folderout 'corr_classic'],8,6,'pdf')
savefig_FC([folderout 'corr_classic'],8,6,'fig')

%% Tracer: find optimal ndt for Correlation function -- dt method
tracklong_tracer = trajs_conc;
for i=1:numel(tracklong_tracer)
    tracklong_tracer(i).X = tracklong_tracer(i).x;
    tracklong_tracer(i).Y = tracklong_tracer(i).y;
    tracklong_tracer(i).Z = tracklong_tracer(i).z;
end

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
yline(0,LineWidth = 2)
%legend('$x$','$y$','$z$','Location','northeast','FontWeight','bold','FontSize',20)
%xlabel('$\tau(s)$','interpreter','latex','FontWeight','bold','FontSize',20)
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$','interpreter','latex','FontWeight','bold','FontSize',20)
grid off;

subplot(2,1,2)
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau.(f),corra.(f)/corra.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end
grid off;
yline(0,LineWidth = 2)
legend('$x$','$y$','$z$','Location','northeast','FontWeight','bold','FontSize',20,'interpreter','latex')
xlabel('$\tau(s)$','interpreter','latex','FontWeight','bold','FontSize',20);
ylabel('$\frac{\langle a(t)a(t+\tau) \rangle} {\langle a^2(t) \rangle}$','interpreter','latex','FontWeight','bold','FontSize',20);
xlim([0 0.04])


folderout = 'corr/';
mkdir(folderout)
savefig_FC([folderout 'corr_dt'],8,6,'pdf')
savefig_FC([folderout 'corr_dt'],8,6,'fig')

save('output_post_processing.mat','Ruu','Raa','Ruufit','Raafit','-append')
end
