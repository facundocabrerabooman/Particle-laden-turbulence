clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/exports/fullg/';

%folderout=folderin;
folderout = '/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/fullg/';
mkdir(folderout)
cd(folderout)

Fs=2996; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Concatenate data
trajs_conc = [];

if 1==pi
    %%% TRACERS


else
    load([folderin filesep 'trajsf_Tracers_01_fullg.mat'],'tracklong')
    trajs_conc = [trajs_conc tracklong];
    1

    % load([folderin filesep 'trajsf_Tracers_02_fullg.mat'],'tracklong')
    % trajs_conc = [trajs_conc tracklong];
    % 
    % load([folderin filesep 'trajsf_Tracers_03_fullg.mat'],'tracklong')
    % trajs_conc = [trajs_conc tracklong];

%%% full g
    % load([folderin filesep 'trajsf_TrCer_1000_13_fullg_tracer.mat'],'tracklong')
    % trajs_conc = [trajs_conc tracklong];
    % 
    % load([folderin filesep 'trajsf_TrCer_1000_14_fullg_tracer.mat'],'tracklong')
    % trajs_conc = [trajs_conc tracklong];
%%%

    load([folderin filesep 'trajsf_Tracers_04_fullg.mat'],'tracklong')
    trajs_conc = [trajs_conc tracklong];
    2

end

Ine=find(arrayfun(@(X)(numel(X.Vx)>4),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);

%%% compute absolute value of vel and acc

for i=1:numel(trajs_conc)
    trajs_conc(i).Vabs = sqrt([trajs_conc(i).Vx].^2+[trajs_conc(i).Vy].^2+[trajs_conc(i).Vz].^2);
    trajs_conc(i).Aabs = sqrt([trajs_conc(i).Ax].^2+[trajs_conc(i).Ay].^2+[trajs_conc(i).Az].^2);
end
%%%

clear tracklong

%%% change reference system
counter=0;
for i=1:numel(trajs_conc)
        trajs_conc(i).X = trajs_conc(i).x;
        trajs_conc(i).Y = trajs_conc(i).z;
        trajs_conc(i).Z = -trajs_conc(i).y;

        trajs_conc(i).Vx = trajs_conc(i).Vx;
        trajs_conc(i).Vy = trajs_conc(i).Vz;
        trajs_conc(i).Vz = -trajs_conc(i).Vy;

        trajs_conc(i).Ax = trajs_conc(i).Ax;
        trajs_conc(i).Ay = trajs_conc(i).Az;
        trajs_conc(i).Az = -trajs_conc(i).Ay;

        %trajs_conc(i).Z = trajs_conc(i).z;
end

%%% If doing microgravity data Only keep data before 5400frames, where the peak in acceleration
%%% happens +
%disp('microgravity?'); pause
% trajs_conc_new=trajs_conc;
% for i=1:numel(trajs_conc)
%     if trajs_conc(i).t(end)>1.8
%         i
%         trajs_conc_new(i)=[];
%     end
% end
% trajs_conc=trajs_conc_new;
disp('GRAVITY goes downwards in z now; the vectors are anti-parallel')


trajs_conc_with_mean_field = trajs_conc;

% try
% save([folderin filesep 'traj_conc'],'trajs_conc','-v7.3')
% catch
% end

clearvars -except folderin folderout fname Fs mycolormap color3 color1 trajs_conc trajs_conc_with_mean_field
%% Get Mean Velocity
if pi==pi
    trajs_conc2=trajs_conc;
    for i=1:numel(trajs_conc)

        trajs_conc2(i).x = trajs_conc(i).x;
        trajs_conc2(i).y = trajs_conc(i).z;
        trajs_conc2(i).z = -trajs_conc(i).y;

        trajs_conc2(i).Yf = trajs_conc(i).Zf;
        trajs_conc2(i).Zf = -trajs_conc(i).Yf;

        trajs_conc2(i).Vy = trajs_conc(i).Vz;
        trajs_conc2(i).Vz = -trajs_conc(i).Vy;
    end
    trajs_conc = trajs_conc2; clear trajs_conc2
    disp('GRAVITY goes downwards in z now; the vectors are anti-parallel')

    dt = [4 6 8 10];
    nbins = [20 20 30];
    threshold = 1;
    gridRange.x = linspace(-30,30,20);
    gridRange.y = linspace(-30,30,20);
    gridRange.z = linspace(-45,35,30);

    [gridsV,meanFieldsV, trajs_conc_minus_mean_field] = meanFields(trajs_conc,Fs,dt,nbins,threshold,1,1,gridRange,1);

    trajs_conc_with_mean_field = trajs_conc;

    %%% Put data minus_mean_field under the same variable names.

    for i=1:numel(trajs_conc_minus_mean_field)
        if ~isnan(trajs_conc_minus_mean_field(i).mVx)
            trajs_conc_minus_mean_field(i).Vx = trajs_conc_minus_mean_field(i).mVx;
            trajs_conc_minus_mean_field(i).Vy = trajs_conc_minus_mean_field(i).mVy;
            trajs_conc_minus_mean_field(i).Vz = trajs_conc_minus_mean_field(i).mVz;
        end
    end

    clear trajs_conc

    %%%

    % try
    % save('trajs_conc_minus_mean_field','trajs_conc_minus_mean_field','-v7.3')
    % catch
    % end
end
%% 1 time - 1 particle statistics
%if pi==pi
%% Calculate & plot velocity and acceleration pdfs

pdfVabs = mkpdf5(trajs_conc_with_mean_field,'Vabs',256,10);
pdfAabs = mkpdf5(trajs_conc_with_mean_field,'Aabs',256,20);

pdfV(1) = mkpdf5(trajs_conc_with_mean_field,'Vx',256,10);
1
pdfV(2) = mkpdf5(trajs_conc_with_mean_field,'Vy',256,10);
2
pdfV(3) = mkpdf5(trajs_conc_with_mean_field,'Vz',256,10);
3

pdfA(1) = mkpdf5(trajs_conc_with_mean_field,'Ax',256,20);
4
pdfA(2) = mkpdf5(trajs_conc_with_mean_field,'Ay',256,20);
5
pdfA(3) = mkpdf5(trajs_conc_with_mean_field,'Az',256,20);
6

% try
% save('output_post_processing.mat','pdfV','pdfA')
% catch end
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

folderout = 'pdfs';
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


folderout = 'pdfs';
mkdir(folderout)
savefig_FC([folderout filesep 'PDF_a'],8,6,'pdf')
savefig_FC([folderout filesep 'PDF_a'],8,6,'fig')

%% Table with moments of distribution
maketable(pdfA,pdfV,pdfVabs,pdfAabs,folderout)

%%






%%%%%%%%%%%%%%%%%% 2 times - 1 particle statistics (Lagrangian statistics)




%% Longitudinal S2

S2L(1)= structFunc_struct(trajs_conc_minus_mean_field,'Vx',2);
S2L(2)= structFunc_struct(trajs_conc_minus_mean_field,'Vy',2);
S2L(3)= structFunc_struct(trajs_conc_minus_mean_field,'Vz',2);

% try
% save('S2L.mat','S2L','-v7.3')
% catch end
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

legend('$S_2^L(x)$','$S_2^L(y)$','$S_2^L(z)$','interpreter','latex',Location='northwest')
ylabel('$S_2^L$','interpreter','latex',FontWeight='bold')
xlabel('$\tau$ (s)','interpreter','latex',FontWeight='bold')
text(8e-4,8e2,'$\tau^2$','interpreter','latex','FontWeight','bold','FontSize',20)
text(1e-2,1.5e4,'$\tau$','interpreter','latex','FontWeight','bold','FontSize',20)
grid on
axis padded

folderout = 'S2L';
mkdir(folderout)
savefig_FC([folderout filesep 'S2L'],8,6,'pdf')
savefig_FC([folderout filesep 'S2L'],8,6,'fig')


%%% add Co inset

axes('Position',[.5 .3 .3 .3])
box on

loglog(S2L(1).tau/Fs,S2L(1).mean./(S2L(1).tau/Fs*5.8e4),'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on
loglog(S2L(2).tau/Fs,S2L(2).mean./(S2L(2).tau/Fs*5.8e4),'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);
loglog(S2L(3).tau/Fs,S2L(3).mean./(S2L(3).tau/Fs*5.8e4),'d',MarkerSize=3,Color=color3(3,:),LineWidth=2);

yline(5.28,'k','LineWidth',3)

ylim([0.08 10])

ylabel('$C_k = S_2^L/(\epsilon \tau)$','interpreter','latex',FontWeight='bold')
xlabel('$\tau$ (s)','interpreter','latex',FontWeight='bold')

savefig_FC([folderout filesep 'S2L_insetC0'],8,6,'pdf')
savefig_FC([folderout filesep 'S2L_insetC0'],8,6,'fig')
%end

%% Eulerian 2-point statistics
%clearvars -except trajs_conc Ine Ine Fs color1 color3 folderin folderout mycolormap

%for j=1:numel(trajs_conc); trajs_conc(j).Tf = trajs_conc(j).t_sec_abs; end % rename Tf field

tic  
[eulerStats, pair] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_minus_mean_field,[0.5 40],40,'off');
toc
%save(['euler' gravity '.mat'],'eulerStats','pair','v7.3')

%clearvars -except eulerStats pair trajs_conc Ine Fs folderin folderout color3 color1
%% Plot
figure;
loglog(eulerStats.r,eulerStats.S2x,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.S2y,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.S2z,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

rS2E = linspace(0.5,40,100);
loglog(rS2E,6e3*rS2E.^1,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(4,4e4,'$r^1$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

folderout = 'S2euler/';
mkdir(folderout)
savefig_custom([folderout 'S2e'],8,6,'pdf')
savefig_custom([folderout 'S2e'],8,6,'fig')

figure;hold on;
semilogy(eulerStats.r,(eulerStats.Splong{2}./2.1).^(3/2)./eulerStats.r'./1e6,'o-')
%semilogy(eulerStats.r,eulerStats.Splong{2},'o-')
grid;
xlabel('r (mm)','Interpreter','latex');
ylabel('$(S_2^E^\parallel / C_2 r^{2/3})^{3/2}$','Interpreter','latex');
set(gca,'FontSize',24);
set(gca,'Xscale','log','Yscale','log');
title('Compensated Eulerian S2ps');
fname = 'compensated_S2_epsilon';
 

folderout = 'S2euler';
mkdir(folderout)
savefig_custom([folderout filesep 'S2el_compensated_epsilon'],8,6,'pdf')
savefig_custom([folderout filesep 'S2el_compensated_epsilon'],8,6,'fig')

%% plot VAt (to check stationary)

figure
t=tiledlayout(4,1,'TileSpacing','tight');
nexttile;
plot(eulerStats.Vmoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.VmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.VmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.VmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
legend('$\langle \sqrt{x^2+y^2+z^2} \rangle$','$\langle x \rangle$','$\langle y \rangle$','$\langle z \rangle$','interpreter','latex',Location='best',FontSize=10);
title('$|V|, \sigma_V, |A|,\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$|V|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Vstd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.VstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.VstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.VstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$\sigma_V$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Amoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.AmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.AmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.AmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$|A|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Astd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.AstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.AstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.AstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/s$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticks(0:Fs/10:size(eulerStats.Astd,2))
xticklabels(num2cell([0:Fs/10:size(eulerStats.Astd,2)]/Fs))
grid on


linkaxes(t.Children,'x')

folderout = 'Vat';
mkdir(folderout)
savefig_custom([folderout filesep 'Vat'],8,6,'pdf')
savefig_custom([folderout filesep 'Vat'],8,6,'fig')

%% plot Spn_abs
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

figure;
loglog(eulerStats.r,eulerStats.SplongAbs{1,1},'d-',MarkerSize=8,Color=color5(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.SplongAbs{1,2},'d-',MarkerSize=8,Color=color5(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,3},'d-',MarkerSize=8,Color=color5(3,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,4},'d-',MarkerSize=8,Color=color5(4,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,5},'d-',MarkerSize=8,Color=color5(5,:),LineWidth=2);

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

folderout = 'SplongAbs';
mkdir(folderout)
savefig_custom([folderout filesep 'SplongAbs'],8,6,'pdf')
savefig_custom([folderout filesep 'SplongAbs'],8,6,'fig')

%% plot Spn
figure;
loglog(eulerStats.r,eulerStats.Splong{1,1},'d-',MarkerSize=8,Color=color5(1,:),LineWidth=1);hold on
loglog(eulerStats.r,eulerStats.Splong{1,2},'d-',MarkerSize=8,Color=color5(2,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,3},'d-',MarkerSize=8,Color=color5(3,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,4},'d-',MarkerSize=8,Color=color5(4,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,5},'d-',MarkerSize=8,Color=color5(5,:),LineWidth=1);

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1,LineWidth=1)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded


folderout = 'Splong';
mkdir(folderout)
savefig_custom([folderout filesep 'Splong'],8,6,'pdf')
savefig_custom([folderout filesep 'Splong'],8,6,'fig')

%% plot Sau
figure
semilogx(eulerStats.r,eulerStats.Sau,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
semilogx(eulerStats.r,eulerStats.Saulong,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=1);hold on

set(gca,FontSize=15)
legend('$S_{au}$','$S_{au}^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

folderout = 'Sau';
mkdir(folderout)
savefig_custom([folderout filesep 'Sau'],8,6,'pdf')
savefig_custom([folderout filesep 'Sau'],8,6,'fig')

%% plot S2E
figure;
loglog(eulerStats.r,eulerStats.S2x,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
loglog(eulerStats.r,eulerStats.S2y,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.S2z,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=1);

loglog(eulerStats.r,7e3*eulerStats.r.^(2/3),'--',Color=color1,LineWidth=1)

set(gca,FontSize=15)
legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(3,3e4,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

folderout = 'S2E';
mkdir(folderout)
savefig_custom([folderout filesep 'S2E'],8,6,'pdf')
savefig_custom([folderout filesep 'S2E'],8,6,'fig')
%% plot epsilon
Ckolomogrov = 2.1;
figure;
loglog(eulerStats.r,(eulerStats.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats.r','d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);

hold on
%figure
loglog(eulerStats.r,abs(eulerStats.Splong{1,3})./(4/5*eulerStats.r)','d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(eulerStats.r,abs(eulerStats.Sau)./2,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

set(gca,FontSize=15)
legend('$(S_2^{\parallel}/C_k)^{3/2}\cdot r^{-1}$','$|S_3^{\parallel}|\cdot (4/5r)^{-1}$','$|S_{au}|/2$','interpreter','latex',Location='best',FontSize=12)
title('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded


folderout = 'epsilon';
mkdir(folderout)
savefig_custom([folderout filesep 'epsilon'],8,6,'pdf')
savefig_custom([folderout filesep 'epsilon'],8,6,'fig')

stop
%% Mean Square Separation
MSD(1) = structFunc_struct(trajs_conc_minus_mean_field,'Xf',2);
1
MSD(2) = structFunc_struct(trajs_conc_minus_mean_field,'Yf',2);
2
MSD(3) = structFunc_struct(trajs_conc_minus_mean_field,'Zf',2);
3

% try
% save('output_post_processing.mat','MSD','-append')
% catch end
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


folderout = 'MSS';
mkdir(folderout)
savefig_FC([folderout filesep  'MSS'],8,6,'pdf')
savefig_FC([folderout filesep 'MSS'],8,6,'fig')

%% Velocity and Acceleration Correlations
if pi==pi
    disp('corr')

    %Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);


    Ruu(1) = xcorr_struct(trajs_conc_minus_mean_field,'Vx',1);
    Ruu(2) = xcorr_struct(trajs_conc_minus_mean_field,'Vy',1);
    Ruu(3) = xcorr_struct(trajs_conc_minus_mean_field,'Vz',1);

    Raa(1) = xcorr_struct(trajs_conc_minus_mean_field,'Ax',1);
    %load('output_post_processing.mat','Ruu','Raa')
    Raa(2) = xcorr_struct(trajs_conc_minus_mean_field,'Ay',1);
    Raa(3) = xcorr_struct(trajs_conc_minus_mean_field,'Az',1);

    %try
    %save('output_post_processing.mat','Ruu','Raa')
    %catch end
    %% Tracer: fit correlation  -- Thomas's
    addpath(genpath('/Users/fcb/Documents/GitHub/Cheng_old'));

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

    folderout = 'corr';
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


    folderout = 'corr';
    mkdir(folderout)
    savefig_FC([folderout 'corr_dt'],8,6,'pdf')
    savefig_FC([folderout 'corr_dt'],8,6,'fig')

    save('output_post_processing.mat','Ruu','Raa','Ruufit','Raafit','-append')
end
