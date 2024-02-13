%% Definitions
clc

addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/'));

folderout = '/Volumes/landau1/TrCer_1000/joint_plots';
cd(folderout)

folder_ddt = '/Volumes/landau1/TrCer_1000/ddt';
folder_fullg = '/Volumes/landau1/TrCer_1000/fullg';
folder_dec = '/Volumes/landau1/TrCer_1000/dec';

%%% Load slip velocities for fullg-ddt-dec

load([folder_ddt filesep 'slipVelCylind_CONC.mat'],'AverSlipVelCylind_conc');
AverSlipVelCylind_conc_ddt = AverSlipVelCylind_conc;

load([folder_fullg filesep 'slipVelCylind_CONC.mat'],'AverSlipVelCylind_conc');
AverSlipVelCylind_conc_fullg = AverSlipVelCylind_conc;

% load([folder_dec filesep 'slipVelCylind_CONC.mat'],'AverSlipVelCylind_conc');
% AverSlipVelCylind_conc_dec = AverSlipVelCylind_conc;

%%% Load particle velocities

% no turbulence case
load('/Volumes/landau1/TrCer_1000/noturb/trajs_TrCer_1000_noturb.mat')
tracklong_noturb = tracklong; clear tracklong

% Fullg particles
load('/Volumes/landau1/TrCer_1000/fullg/tracklongP_conc.mat')
tracklong_fullg = tracklongP_conc; clear tracklongP_conc

% DDT particles
load('/Volumes/landau1/TrCer_1000/ddt/tracklongP_conc.mat')
tracklong_ddt = tracklongP_conc; clear tracklongP_conc

% % Dec particles
% load('/Volumes/landau1/TrCer_1000/dec/tracklongP_conc.mat')
% tracklong_dec = tracklongP_conc; clear tracklongP_conc

%% PDFs

%load([folderout filesep 'tracklongP_conc'],'tracklongP_conc')
%tracklong=tracklongP_conc; clear tracklongP_conc

%pdfV(2) = mkpdf5(tracklong,'Vy',256,10);
%pdfV(2) = mkpdf5(tracklong,'Vy',100,5);

pdfSV(1) = mkpdf5(AverSlipVelCylind_conc_fullg,'Urelmean_vert',100,10);
1
pdfSV(2) = mkpdf5(AverSlipVelCylind_conc_ddt,'Urelmean_vert',100,10);
2
%pdfSV(3) = mkpdf5(AverSlipVelCylind_conc_dec,'Urelmean_vert',100,5);
%3

pdfV(1) = mkpdf5(tracklong_fullg,'Vy',100,10);
pdfV(2) = mkpdf5(tracklong_ddt,'Vy',100,10);
%pdfV(3) = mkpdf5(tracklong_dec,'Vy',100,10);

%% Plot PDFs

figure(1); hold on; clf

%%% plot slip vel
semilogy(pdfSV(1).xpdf,pdfSV(1).pdf,'r-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfSV(1).mean,'r',LineWidth=3)

semilogy(pdfSV(2).xpdf,pdfSV(2).pdf,'g-o',MarkerSize=5,LineWidth=2);
xline(pdfSV(2).mean,'g',LineWidth=3)

% semilogy(pdfSV(3).xpdf,pdfSV(3).pdf,'b-o',MarkerSize=5,LineWidth=2);
% xline(pdfSV(3).mean,'b',LineWidth=3)

%%% plot particle vel
semilogy(pdfV(1).xpdf,pdfV(1).pdf,'m-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfV(1).mean,'m',LineWidth=3)

semilogy(pdfV(2).xpdf,pdfV(2).pdf,'c-o',MarkerSize=5,LineWidth=2);
xline(pdfV(2).mean,'c',LineWidth=3)

% semilogy(pdfV(3).xpdf,pdfV(3).pdf,'y-o',MarkerSize=5,LineWidth=2);
% xline(pdfV(3).mean,'y',LineWidth=3)

%%% plot no turb vel
xline(mean(vertcat(tracklong_noturb(:).Vy)),'k',LineWidth=3)


%legend('$SlipV-FULLG$','$SlipV-DDT$','$SlipV-DEC$','FULLG','DDT','DEC','No Turbulence - full g','interpreter','latex',Location='northeast');
legend('$SlipV-FULLG$','Mean','$SlipV-DDT$','Mean','FULLG','Mean','DDT','Mean','No Turbulence - full g','interpreter','latex',Location='northeast');
%ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
xlabel('Vertical Velocity (mm/s)',FontWeight='bold')
xticks(-800:100:800)
grid on
box on
axis padded


folderout_tmp = 'pdfs';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'PDF_v'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'PDF_v'],8,6,'fig')

%% Plot mean velocities

figure(2);hold on; clf

%%% fullg-ddt-dec slip settling velocity
mean2 = pdfSV(1).mean;
std = pdfSV(1).std/2;
yline(mean2,'r',LineWidth=3)
% vertices = [0, mean2-std; 1, mean2-std; 1, mean2+std; 0, mean2+std];
% patch('Vertices', vertices, 'Faces', [1 2 3 4],'EdgeColor','none', 'FaceColor', 'r', 'FaceAlpha', 0.2, 'LineWidth', 2);

mean2 = pdfSV(2).mean;
std = pdfSV(2).std/2;
yline(mean2,'g',LineWidth=3)
% vertices = [0, mean2-std; 1, mean2-std; 1, mean2+std; 0, mean2+std];
% patch('Vertices', vertices, 'Faces', [1 2 3 4],'EdgeColor','none', 'FaceColor', 'g', 'FaceAlpha', 0.2, 'LineWidth', 2);

% mean = pdfSV(3).mean2;
% std = pdfSV(3).std2/2;
% yline(mean,'b',LineWidth=3)
% vertices = [0, mean2-std; 1, mean2-std; 1, mean2+std; 0, mean2+std];
% patch('Vertices', vertices, 'Faces', [1 2 3 4],'EdgeColor','none', 'FaceColor', 'b', 'FaceAlpha', 0.2, 'LineWidth', 2);

%%% fullg-ddt-dec particle settling velocity
mean2 = pdfV(1).mean;
std = pdfV(1).std/2;
yline(mean2,'m',LineWidth=3)
% vertices = [0, mean2-std; 1, mean2-std; 1, mean2+std; 0, mean2+std];
% patch('Vertices', vertices, 'Faces', [1 2 3 4],'EdgeColor','none', 'FaceColor', 'm', 'FaceAlpha', 0.2, 'LineWidth', 2);

mean2= pdfV(2).mean;
std = pdfV(2).std/2;
yline(mean2,'c',LineWidth=3)
% vertices = [0, mean2-std; 1, mean2-std; 1, mean2+std; 0, mean2+std];
% patch('Vertices', vertices, 'Faces', [1 2 3 4],'EdgeColor','none', 'FaceColor', 'c', 'FaceAlpha', 0.2, 'LineWidth', 2);

% mean2 = pdfV(3).mean;
% std = pdfV(3).std/2;
% yline(mean2,'y',LineWidth=3)
% vertices = [0, mean2-std; 1, mean2-std; 1, mean2+std; 0, mean2+std];
% patch('Vertices', vertices, 'Faces', [1 2 3 4],'EdgeColor','none', 'FaceColor', 'y', 'FaceAlpha', 0.2, 'LineWidth', 2);

%%% plot no turb vel
yline(mean(vertcat(tracklong_noturb.Vy)),'k',LineWidth=3)


%legend('$SlipV-FULLG$','$SlipV-DDT$','$SlipV-DEC$','FULLG','DDT','DEC','No Turbulence - full g','interpreter','latex',Location='northeast');
%ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
ylabel('Vertical Velocity (mm/s)',FontWeight='bold')
box on; grid on

folderout_tmp = 'meanvels';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'meanvel'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'meanvel'],8,6,'fig')

%% Compare PDF shapes - plot normlized distributions

figure(3); hold on; clf

%%% plot slip vel
semilogy(pdfSV(1).xpdfn,pdfSV(1).pdfn,'r-o',MarkerSize=5,LineWidth=2);hold on

semilogy(pdfSV(2).xpdfn,pdfSV(2).pdfn,'g-o',MarkerSize=5,LineWidth=2);

% semilogy(pdfSV(3).xpdfn,pdfSV(3).pdfn,'b-o',MarkerSize=5,LineWidth=2);

%%% plot particle vel
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'m-o',MarkerSize=5,LineWidth=2);hold on

semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'c-o',MarkerSize=5,LineWidth=2);

%semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'y-o',MarkerSize=5,LineWidth=2);


%legend('$SlipV-FULLG$','$SlipV-DDT$','$SlipV-DEC$','FULLG','DDT','DEC','No Turbulence - full g','interpreter','latex',Location='northeast');
legend('$SlipV-FULLG$','$SlipV-DDT$','FULLG','DDT','interpreter','latex',Location='south');
ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
xlabel('Vertical Velocity (mm/s)',FontWeight='bold')
xticks(-10:2:10)
grid on
box on
axis padded


savefig_FC([folderout_tmp filesep 'PDF_comparison'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'PDF_comparison'],8,6,'fig')
%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BOTH: slip velocity averaged over experiments

% Initialize variables for storing sum and mean of relative velocities
nbinr = 30;
nbint = 30;
nidxmax = 1000;

for kexp = 1:Nexp
    Urel_allMatrices.x{:,kexp} = zeros(nbinr,nbint,nidxmax)*NaN;
    Urel_allMatrices.y{:,kexp} = zeros(nbinr,nbint,nidxmax)*NaN;
    Urel_allMatrices.z{:,kexp} = zeros(nbinr,nbint,nidxmax)*NaN;
end
fieldname = fieldnames(Urel_allMatrices);

% Create a spherical mesh
[rr,tt,~,~,X,Y] = SphericalMesh([Rmin Rmax],nbinr,nbint);

% Loop over different experiments (kexp)
for kexp = 1:Nexp
    % Load slip velocity data for the current experiment
    fname = ['slipVeloCylind_',num2str(kexp),'.mat'];
    load([fpathP_slipVeloData fname],'slipVeloCylind')

    % Extract relevant columns from slip velocity data
    rho = vertcat(slipVeloCylind.rho);
    theta = vertcat(slipVeloCylind.theta);
    z = vertcat(slipVeloCylind.z);
    urel = vertcat(slipVeloCylind.Urel);

    Urel.x = urel(:,1);
    Urel.y = urel(:,2);
    Urel.z = urel(:,3);
    %     Urel.norm = cell2mat(arrayfun(@(X) sqrt(X.x.^2+X.y.^2+X.z.^2),Urel,'UniformOutput',false));

    % Transform cylindrical coordinates to 2D polar coordinates
    rho2D = sqrt(rho.^2 + z.^2);
    theta2D = atan2(rho,z);

    % Iterate over each bin and calculate the sum of Urel
    for i = 1:nbinr-1
        for j = 1:nbint-1
            ind = find(rho2D >= rr(i) & rho2D < rr(i+1) & theta2D >= tt(j) & theta2D< tt(j+1));
            if numel(ind)~=0
                for kf = 1:numel(fieldname)
                    extendedVector = nan(1, nidxmax);
                    extendedVector(1:numel(ind)) = Urel.(fieldname{kf})(ind);
                    Urel_allMatrices.(fieldname{kf}){:,kexp}(j,i,:) = extendedVector;
                end
            end

        end
    end
    clear slipVeloCylind
end

for i = 1:numel(fieldname)
    [averMeanUrel.(fieldname{i}), averStdUrel.(fieldname{i}),countUrel] = averUrelMap(Urel_allMatrices.(fieldname{i}));
end
clear Urel_allMatrices

save([fpathP_slipVeloData 'averUrel.mat'],'averMeanUrel','averStdUrel','countUrel')

%% Plot the binned mean values using pseudocolor

etaK=1; % FCB

etaKMMS = etaK*1e3; % mm

etaKMMS = 1;
Vpg_tsABS = 1; % absValue of the case without turburlence, here we use mm/s

subscripts = {'\rho','\theta','z'};
fnamesub = {'rho','theta','z'};
for i = 1:2
    for kf = 1:numel(fieldname)
        if i == 1
            averUrel = averMeanUrel;
            colstr = ['$\langle U^{rel}_{' subscripts{kf} '}\rangle/\langle V^{Re=0}_{g} \rangle$'];
            fnamestr = 'mean';
            caxisLim = [-1 0];
        else
            averUrel = averStdUrel;
            colstr = ['$\sigma( U^{rel}_{' subscripts{kf} '})/\langle V^{Re=0}_{g} \rangle$'];
            fnamestr = 'std';
            caxisLim = [0 1];
        end
        figure;
        pcolor(X/etaKMMS,Y/etaKMMS, averUrel.(fieldname{kf})/Vpg_tsABS);
        xlabel('$z/\eta_K$');
        ylabel('$r/\eta_K$');
        %     colormap(parula(32))
        col =colorbar;
        %caxis(caxisLim)
        hold on
        rectangle('Position', [[0, 0] - 0.5, 2*0.5, 2*0.5]/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'b','FaceColor','k');
        quiver(0, 0, 4/etaKMMS, 0, 0, 'r', 'LineWidth', 6);
        %     rectangle('Position', [[0, 0] - 1, 2*1, 2*1]*2/etaK, 'Curvature', [1, 1], 'EdgeColor', 'r','LineStyle','--',LineWidth=2);
        ylabel(col,colstr,'interpreter','latex')
        col.TickLabelInterpreter = "latex";
        %  savefig_FC([fpathP_slipVeloData filesep fnamestr 'Urel_' fnamesub{kf}],8,6,'pdf')
        stop
    end
end

figure;
pcolor(X/etaKMMS,Y/etaKMMS, countUrel);
xlabel('$z/\eta_K$');
ylabel('$r/\eta_K$');
% colormap(parula(32))
col =colorbar;
caxis([100 5000])
hold on
rectangle('Position', [[0, 0] - 0.5, 2*0.5, 2*0.5]/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'b','FaceColor','k');
quiver(0, 0, 4/etaKMMS, 0, 0, 'r', 'LineWidth', 6);
rectangle('Position', [[0, 0] - 1, 2*1, 2*1]*2/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'r','LineStyle','--',LineWidth=2);
ylabel(col,'$N_{count}$','interpreter','latex')
col.TickLabelInterpreter = "latex";

savefig_FC([fpathP_slipVeloData filesep 'Urel_count'],8,6,'pdf')
savefig_FC([fpathP_slipVeloData filesep 'Urel_count'],8,6,'fig')
