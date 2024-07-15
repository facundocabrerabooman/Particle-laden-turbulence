%% Definitions
clfc

addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/'));

folderout = '/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Manuscripts/rsi_ddt/inertial_cer/g&ug';
mkdir(folderout)
cd(folderout)

%folder_ddt = '/Volumes/landau1/TrCer_analysis_paper#1/exports/ddt_pairs/';
%folder_fullg = '/Volumes/landau1/TrCer_analysis_paper#1/exports/fullg_pairs/';
%folder_dec = '';

%%% Load slip velocities for fullg-ddt-dec

load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/slipVeloData/slipVeloData_R_10/slipVelCylind_ddt_CONC.mat','AverSlipVelCylind_conc');
AverSlipVelCylind_conc_ddt_R10 = AverSlipVelCylind_conc; clear AverSlipVelCylind_conc

load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/slipVeloData/slipVeloData_R_10/slipVelCylind_fullg_CONC.mat','AverSlipVelCylind_conc');
AverSlipVelCylind_conc_fullg_R10 = AverSlipVelCylind_conc; clear AverSlipVelCylind_conc


%%% Load particle velocities

% no turbulence case
load('/Volumes/landau1/TrCer_analysis_paper#1/noturb/trajs_TrCer_1000_noturb.mat')
tracklong_noturb = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong_noturb)==1);
tracklong_noturb = tracklong_noturb(Ine);

% Fullg particles
load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/particle_fullg_conc.mat')
%load('/Users/fcb/Downloads/particles/fullg/trajsf_TrCer_1000_09_fullg_particle.mat')
tracklong_fullg = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong_fullg)==1);
tracklong_fullg = tracklong_fullg(Ine);

% DDT particles
load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/particle_ddt_conc.mat')
tracklong_ddt = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong_ddt)==1);
tracklong_ddt = tracklong_ddt(Ine);

% Dec particles
% load('/Volumes/landau1/TrCer_1000/dec/tracklongP_conc.mat')
% tracklong_dec = tracklong; clear tracklong


%% Substract Mean Flow
clc
load('/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/fullg/MeanFlow_Interpolant.mat')
Fx_fullg = Fx;
Fy_fullg = Fy;
Fz_fullg = Fz;

for i=1:numel(tracklong_fullg)

    MeanVel_x = Fx_fullg(tracklong_fullg(i).Xf,tracklong_fullg(i).Yf,tracklong_fullg(i).Zf);
    MeanVel_y = Fy_fullg(tracklong_fullg(i).Xf,tracklong_fullg(i).Yf,tracklong_fullg(i).Zf);
    MeanVel_z = Fz_fullg(tracklong_fullg(i).Xf,tracklong_fullg(i).Yf,tracklong_fullg(i).Zf);

    tracklong_fullg(i).Vx_minus_meanflow = tracklong_fullg(i).Vx - MeanVel_x;
    tracklong_fullg(i).Vy_minus_meanflow = tracklong_fullg(i).Vy - MeanVel_y;
    tracklong_fullg(i).Vz_minus_meanflow = tracklong_fullg(i).Vz - MeanVel_z;


    clear MeanVel_x MeanVel_y MeanVel_z
end


load('/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/ddt/MeanFlow_Interpolant.mat')
Fx_ddt = Fx;
Fy_ddt = Fy;
Fz_ddt = Fz; clear Fx Fy Fz

for i=1:numel(tracklong_ddt)

    MeanVel_x = Fx_ddt(tracklong_ddt(i).Xf,tracklong_ddt(i).Yf,tracklong_ddt(i).Zf);
    MeanVel_y = Fy_ddt(tracklong_ddt(i).Xf,tracklong_ddt(i).Yf,tracklong_ddt(i).Zf);
    MeanVel_z = Fz_ddt(tracklong_ddt(i).Xf,tracklong_ddt(i).Yf,tracklong_ddt(i).Zf);

    tracklong_ddt(i).Vx_minus_meanflow = tracklong_ddt(i).Vx - MeanVel_x;
    tracklong_ddt(i).Vy_minus_meanflow = tracklong_ddt(i).Vy - MeanVel_y;
    tracklong_ddt(i).Vz_minus_meanflow = tracklong_ddt(i).Vz - MeanVel_z;

    clear MeanVel_x MeanVel_y MeanVel_z
end

%%%%%% DEC %%%%%%%%
if 1==pi
    %load('/Volumes/landau1/Tracers/dec_13-14-15-16/MeanFlow_Interpolant.mat')
    load('/Volumes/landau1/Tracers/fullg_13-14-15-16/MeanFlow_Interpolant.mat') % I do not have data for tracers in dec so I use fullg 3-25-24
    Fx_ddt = Fx;
    Fy_ddt = Fy;
    Fz_ddt = Fz; clear Fx Fy Fz

    for i=1:numel(tracklong_dec)

        MeanVel_x = Fx_ddt(tracklong_dec(i).Xf,tracklong_dec(i).Yf,tracklong_dec(i).Zf);
        MeanVel_y = Fy_ddt(tracklong_dec(i).Xf,tracklong_dec(i).Yf,tracklong_dec(i).Zf);
        MeanVel_z = Fz_ddt(tracklong_dec(i).Xf,tracklong_dec(i).Yf,tracklong_dec(i).Zf);

        tracklong_dec(i).Vx_minus_meanflow = tracklong_dec(i).Vx - MeanVel_x;
        tracklong_dec(i).Vy_minus_meanflow = tracklong_dec(i).Vy - MeanVel_y;
        tracklong_dec(i).Vz_minus_meanflow = tracklong_dec(i).Vz - MeanVel_z;

        clear MeanVel_x MeanVel_y MeanVel_z
    end
end
% get rid of NaN mean vels

tracklong_ddt = tracklong_ddt(~arrayfun(@(x) any(isnan(x.Vy_minus_meanflow)), tracklong_ddt));

tracklong_fullg = tracklong_fullg(~arrayfun(@(x) any(isnan(x.Vy_minus_meanflow)), tracklong_fullg));

%tracklong_dec = tracklong_dec(~arrayfun(@(x) any(isnan(x.Vy_minus_meanflow)), tracklong_dec));

% compute distributions
pdf_Vvert_minus_meanflow(1) = mkpdf5(tracklong_fullg,'Vy_minus_meanflow',100,10);
pdf_Vvert_minus_meanflow(2) = mkpdf5(tracklong_ddt,'Vy_minus_meanflow',100,10);
%pdf_Vvert_minus_meanflow(3) = mkpdf5(tracklong_dec,'Vy_minus_meanflow',100,10);
pdfV(1) = mkpdf5(tracklong_fullg,'Vy',100,10);
pdfV(2) = mkpdf5(tracklong_ddt,'Vy',100,10);
%pdfV(3) = mkpdf5(tracklong_dec,'Vy',100,10);

%% Plot Slip Velocity Signal

figure(10);hold on; box on
vel3d=[];
start=1;
finish=100;

counter=0;
for i=start:finish
    counter=counter+1;
    vel3d(counter,:) = AverSlipVelCylind_conc_ddt_R10(i).Urelmean;
end

plot(vertcat(AverSlipVelCylind_conc_ddt_R10(start:finish).t)./2996,vel3d(:,1),'r.-')
%plot(vertcat(AverSlipVelCylind_conc_ddt_R10(start:finish).t),vel3d(:,1),'r.-')
plot(vertcat(AverSlipVelCylind_conc_ddt_R10(start:finish).t)./2996,vel3d(:,2),'g.-')
plot(vertcat(AverSlipVelCylind_conc_ddt_R10(start:finish).t)./2996,vel3d(:,3),'b.-')

legend({'x','y (g)','z'})
xlabel('t (s)')
ylabel('Vel (mm/s)')
grid on;box on

savefig_FC('velocity_signal',8,6,'pdf')
savefig_FC('velocity_signal',8,6,'fig')

%% Compute fluctuations

vel3d=[];

for i=1:numel(AverSlipVelCylind_conc_ddt_R10)
    vel3d(i,:) = AverSlipVelCylind_conc_ddt_R10(i).Urelmean;
end

sigmax = var(vel3d(:,1))
sigmay = var(vel3d(:,2))
sigmaz = var(vel3d(:,3))

%%% fullg

vel3d_fullg=[];

for i=1:numel(AverSlipVelCylind_conc_fullg_R10)
    vel3d_fullg(i,:) = AverSlipVelCylind_conc_fullg_R10(i).Urelmean;
end

sigmax = var(vel3d_fullg(:,1))
sigmay = var(vel3d_fullg(:,2))
sigmaz = var(vel3d_fullg(:,3))
%% Vertical velocity versus time bin

figure(1); clf; hold on; grid on; box on

%%% Fullg
tbins_fullg = (0:.5:4);
vels_in_each_bin = cell(1, length(tbins_fullg)-1);
meanvel_in_each_bin = zeros(1, length(tbins_fullg)-1);

meanVy = cell(1, length(tbins_fullg)-1);
countVy = cell(1, length(tbins_fullg)-1);

Vtbins_fullg = cell(1, length(tbins_fullg)-1);

for i = 1:numel(tracklong_fullg)
    %%%%% Quantity to plot
    Yf = tracklong_fullg(i).Tf;
    %  Vy = tracklong_fullg(i).Vy_minus_meanflow;
    Vy = tracklong_fullg(i).Vy;
    %%%%%

    [~, binIdx] = histc(Yf, tbins_fullg);

    binIdx(binIdx == 0) = 1;
    binIdx(binIdx == length(tbins_fullg)) = length(tbins_fullg) - 1;

    for j = 1:length(tbins_fullg)-1
        %%%% Decide if removing mean vel or not.
        Vtbins_fullg{j} = [Vtbins_fullg{j}; Vy(binIdx == j)];
        %%%%
    end
end

for j = 1:length(tbins_fullg)-1
    meanVy{j} = mean(Vtbins_fullg{j});
    countVy{j} = numel(Vtbins_fullg{j});
end

meanVy_fullg = cell2mat(meanVy);
countVy_fullg = cell2mat(countVy);


%%% DDT
tbins_ddt = (0:0.5:4);
vels_in_each_bin = cell(1, length(tbins_ddt)-1);
meanvel_in_each_bin = zeros(1, length(tbins_ddt)-1);

meanVy = cell(1, length(tbins_ddt)-1);
countVy = cell(1, length(tbins_ddt)-1);

Vtbins_ddt = cell(1, length(tbins_ddt)-1);

for i = 1:numel(tracklong_ddt)
    %%%%% Quantity to plot
    Yf = tracklong_ddt(i).Tf;
    % Vy = tracklong_ddt(i).Vy_minus_meanflow;
    Vy = tracklong_ddt(i).Vy;
    %%%%%

    [~, binIdx] = histc(Yf, tbins_ddt);

    binIdx(binIdx == 0) = 1;
    binIdx(binIdx == length(tbins_ddt)) = length(tbins_ddt) - 1;

    for j = 1:length(tbins_ddt)-1
        %%%% Decide if removing mean vel or not.
        Vtbins_ddt{j} = [Vtbins_ddt{j}; Vy(binIdx == j)];
        %%%%
    end
end

for j = 1:length(tbins_ddt)-1
    meanVy{j} = mean(Vtbins_ddt{j});
    countVy{j} = numel(Vtbins_ddt{j});
end

meanVy_ddt = cell2mat(meanVy);
countVy_ddt = cell2mat(countVy);

if 1==pi
    %%% DEC
    tbins_dec = (0:0.5:4);
    vels_in_each_bin = cell(1, length(tbins_dec)-1);
    meanvel_in_each_bin = zeros(1, length(tbins_dec)-1);

    meanVy = cell(1, length(tbins_dec)-1);
    countVy = cell(1, length(tbins_dec)-1);

    Vtbins_dec = cell(1, length(tbins_dec)-1);

    for i = 1:numel(tracklong_dec)
        %%%%% Quantity to plot
        Yf = tracklong_dec(i).Tf;
        Vy = tracklong_dec(i).Vy_minus_meanflow;
        %Vy = tracklong_dec(i).Vy;
        %%%%%

        [~, binIdx] = histc(Yf, tbins_dec);

        binIdx(binIdx == 0) = 1;
        binIdx(binIdx == length(tbins_dec)) = length(tbins_dec) - 1;

        for j = 1:length(tbins_dec)-1
            %%%% Decide if removing mean vel or not.
            Vtbins_dec{j} = [Vtbins_dec{j}; Vy(binIdx == j)];
            %%%%
        end
    end

    for j = 1:length(tbins_dec)-1
        meanVy{j} = mean(Vtbins_dec{j});
        countVy{j} = numel(Vtbins_dec{j});
    end

    meanVy_dec = cell2mat(meanVy);
    countVy_dec = cell2mat(countVy);
end

% plot
video = VideoWriter('vel_vs_time','MPEG-4');
video.FrameRate = 3e-1;
open(video);

savefig_FC('vertical_vel_versus_time',8,6,'pdf')
xlabel('t (sec)');
ylabel('Mean Vert Vel (mm/s)');
title('Mean Value of Vert. Vel. minus Mean Flow');
xlim([-0.5 5])

% Initialize the loop
for i = 1:length(tbins_fullg)-1
    % Update the plot with new data

    bar(tbins_fullg(i), meanVy_fullg(i),0.4,'r','FaceAlpha', 0.8,'EdgeColor', 'k', 'LineWidth', 1);
    hold on;

    % Pause for a short duration to visualize the progression
    pause(.1);

    % Use drawnow to update the plot immediately
    drawnow;

    % Write the current frame to the video
    frame = getframe(gcf);
    writeVideo(video, frame);

    % Duplicate frames to simulate a lower frame rate
    % for j = 1:10
    %     % Write the current frame to the video
    %     frame = getframe(gcf);
    %     writeVideo(video, frame);
    % end
end

for i = 1:length(tbins_ddt)-1
    % Update the plot with new data

    nan_indices1 = find(isnan(meanVy_fullg));
    bar(tbins_ddt(i)+tbins_fullg(nan_indices1(1)), meanVy_ddt(i),0.4,'g','FaceAlpha', 0.8,'EdgeColor', 'k', 'LineWidth', 1);
    hold on;

    % Pause for a short duration to visualize the progression
    pause(0.1);

    % Use drawnow to update the plot immediately
    drawnow;

    % Write the current frame to the video
    frame = getframe(gcf);
    writeVideo(video, frame);
end

%%%%%%%%%% DEC %%%%%%%%%%
if 1==pi
    for i = 1:length(tbins_dec)-1
        % Update the plot with new data
    
        nan_indices2 = find(isnan(meanVy_ddt));
        bar(tbins_dec(i)+tbins_fullg(nan_indices1(1))+tbins_ddt(nan_indices2(1)), meanVy_dec(i),0.4,'b','FaceAlpha', 0.8,'EdgeColor', 'k', 'LineWidth', 1);
        % Pause for a short duration to visualize the progression
    
        pause(0.1);
    
        % Use drawnow to update the plot immediately
        drawnow;
    
        % Write the current frame to the video
        frame = getframe(gcf);
        writeVideo(video, frame);
    end
end
close(video)

savefig_FC('vertical_vel_versus_time',8,6,'pdf')
savefig_FC('vertical_vel_versus_time',8,6,'fig')

%% Plot
figure(10);clf

%%% fullg-ddt-dec particle settling velocity
mean2 = pdfV(1).mean;
%mean2 = mean(vertcat(tracklong_fullg.Vy_minus_meanflow));
h1=yline(mean2,'m',LineWidth=3);

mean2= pdfV(2).mean;
h2=yline(mean2,'c',LineWidth=3);

%%% fullg-ddt-dec particle settling velocity minus mean flow
mean2 = pdf_Vvert_minus_meanflow(1).mean;
h3=yline(mean2,'-.m',LineWidth=3);

mean2= pdf_Vvert_minus_meanflow(2).mean;
h4=yline(mean2,'-.c',LineWidth=3);

%%% plot no turb vel
h5=yline(mean(vertcat(tracklong_noturb.Vy)),'k',LineWidth=3);

ylabel('Vertical Velocity (mm/s)',FontWeight='bold')
box on; grid on
legend([h1 h2 h3 h4 h5],'FULLG','DDT','FULLG MINUS MEANFIELD','DDT MINUS MEANFIELD','No Turbulence - full g','interpreter','latex',Location='northeast');


folderout_tmp = 'vels_raw_minusmeanfield';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'meanvel'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'meanvel'],8,6,'fig')

%% Regions that particles explore.

for i=1:numel(tracklong_fullg)

    pos_tmp = [tracklong_fullg(i).Xf,tracklong_fullg(i).Yf,tracklong_fullg(i).Zf];
    tracklong_fullg(i).r = sqrt(sum(pos_tmp.^2, 2));

end


for i=1:numel(tracklong_ddt)

    pos_tmp = [tracklong_ddt(i).Xf,tracklong_ddt(i).Yf,tracklong_ddt(i).Zf];
    tracklong_ddt(i).r = sqrt(sum(pos_tmp.^2, 2));

end


%%%%%%%%%%%%%%%% plot
figure(11);clf

scatter3(vertcat(tracklong_fullg.Xf),vertcat(tracklong_fullg.Yf),vertcat(tracklong_fullg.Zf),10,'bo','filled'); hold on
scatter3(vertcat(tracklong_ddt.Xf),vertcat(tracklong_ddt.Yf),vertcat(tracklong_ddt.Zf),10,'mo','filled')

ylabel('Vertical (mm)',FontWeight='bold')
zlabel('Z (mm)',FontWeight='bold')
xlabel('X (mm)',FontWeight='bold')
box on; grid on
legend('FULLG','DDT','interpreter','latex',Location='northeast');

folderout_tmp = 'traj3D';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'traj3D'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'traj3D'],8,6,'fig')



%% PDFs

%%% Create variable with velocity to use mkpdf5
for j=1:numel(AverSlipVelCylind_conc_fullg_R10)
    AverSlipVelCylind_conc_fullg_R10(j).VerticalVel = AverSlipVelCylind_conc_fullg_R10(j).Urel(:,2);
end
for j=1:numel(AverSlipVelCylind_conc_ddt_R10)
    AverSlipVelCylind_conc_ddt_R10(j).VerticalVel = AverSlipVelCylind_conc_ddt_R10(j).Urel(:,2);
end


pdfSV(1) = mkpdf5(AverSlipVelCylind_conc_fullg_R10,'VerticalVel',100,10);
1
pdfSV(2) = mkpdf5(AverSlipVelCylind_conc_ddt_R10,'VerticalVel',100,10);
2
%pdfSV(3) = mkpdf5(AverSlipVelCylind_conc_dec,'Urelmean_vert',100,5);
%3

pdfV(1) = mkpdf5(tracklong_fullg,'Vy',100,10);
pdfV(2) = mkpdf5(tracklong_ddt,'Vy',100,10);
%pdfV(3) = mkpdf5(tracklong_dec,'Vy',100,10);




%% Plot PDFs

figure(1); hold on; clf

if pi==pi %dont have that data yet

    %%% plot slip vel
    semilogy(pdfSV(1).xpdf,pdfSV(1).pdf,'r-o',MarkerSize=5,LineWidth=2);hold on
    xline(pdfSV(1).mean,'r',LineWidth=3)

    semilogy(pdfSV(2).xpdf,pdfSV(2).pdf,'g-o',MarkerSize=5,LineWidth=2);hold on
    xline(pdfSV(2).mean,'g',LineWidth=3)

    % semilogy(pdfSV(3).xpdf,pdfSV(3).pdf,'b-o',MarkerSize=5,LineWidth=2);
    % xline(pdfSV(3).mean,'b',LineWidth=3)
end

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

if pi==pi %dont have that data yet

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
end
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
legend('$SlipV-FULLG$','$SlipV-DDT$','FULLG','DDT','No Turbulence - full g','interpreter','latex',Location='northeast');
ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
ylabel('Vertical Velocity (mm/s)',FontWeight='bold')
box on; grid on


folderout_tmp = 'meanvels';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'meanvel'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'meanvel'],8,6,'fig')

%% Compare PDF shapes - plot normalized distributions

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

%%
figure(4);clf;clc

Rtmp = [2 4 6 10 15:10:75];
for i=1:11
    i/11

    try
        % fullg
        load(['/Volumes/landau1/TrCer_1000/fullg/slipVeloData/slipVeloData_R_' num2str(Rtmp(i)) '/slipVelCylind_CONC.mat'],'AverSlipVelCylind_conc');
        for j=1:numel(AverSlipVelCylind_conc) % add new field to use mkfpdf5
            AverSlipVelCylind_conc(j).VerticalVel = AverSlipVelCylind_conc(j).Urel(:,2);
        end

        pdfSV_tmp = mkpdf5(AverSlipVelCylind_conc,'VerticalVel',100,10);
        scatter(Rtmp(i),pdfSV_tmp.mean,200,'ro','filled'); hold on
        % ddt
        load(['/Volumes/landau1/TrCer_1000/ddt/slipVeloData/slipVeloData_R_' num2str(Rtmp(i)) '/slipVelCylind_CONC.mat'],'AverSlipVelCylind_conc');
        for j=1:numel(AverSlipVelCylind_conc) % add new field to use mkfpdf5
            AverSlipVelCylind_conc(j).VerticalVel = AverSlipVelCylind_conc(j).Urel(:,2);
        end

        pdfSV_tmp = mkpdf5(AverSlipVelCylind_conc,'VerticalVel',100,10);
        scatter(Rtmp(i),pdfSV_tmp.mean,200,'gpentagram','filled'); hold on
    catch end

end

% Particles' raw velocity
pdfV(1) = mkpdf5(tracklong_fullg,'Vy',100,10);
pdfV(2) = mkpdf5(tracklong_ddt,'Vy',100,10);
% fullg
mean2 = pdfV(1).mean;
yline(mean2,'m',LineWidth=3)
% ddt
mean2= pdfV(2).mean;
yline(mean2,'c',LineWidth=3)

% no turb
yline(mean(vertcat(tracklong_noturb.Vy)),'k',LineWidth=3)

box on
xlabel('R')
ylabel('Vertical Velocity (mm/s)',FontWeight='bold')
grid on


folderout_tmp = 'Different_Rs';
mkdir(folderout_tmp)

savefig_FC([folderout_tmp filesep 'vel_vs_R'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'vel_vs_R'],8,6,'fig')
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

fpathP_slipVeloData = '/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/slipVeloData/slipVeloData_R_10/';

Rmin = 0;
Rmax = 10;

kexp=1;
Nexp = 1;

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
    fname = 'slipVelCylind_fullg_CONC.mat';
    load([fpathP_slipVeloData fname],'AverSlipVelCylind_conc')

    % Extract relevant columns from slip velocity data
    rho = vertcat(AverSlipVelCylind_conc_fullg_R10.rho);
    theta = vertcat(AverSlipVelCylind_conc_fullg_R10.theta);
    z = vertcat(AverSlipVelCylind_conc_fullg_R10.z);
    urel = vertcat(AverSlipVelCylind_conc_fullg_R10.Urel);

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

save('averUrel.mat','averMeanUrel','averStdUrel','countUrel')

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
         savefig_FC([folderout filesep fnamestr 'Urel_' fnamesub{kf}],8,6,'pdf')
        
    end
end

figure;
pcolor(X/etaKMMS,Y/etaKMMS, countUrel);
xlabel('$z/\eta_K$');
ylabel('$r/\eta_K$');
% colormap(parula(32))
col =colorbar;
caxis([0 500])
hold on
rectangle('Position', [[0, 0] - 0.5, 2*0.5, 2*0.5]/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'b','FaceColor','k');
quiver(0, 0, 4/etaKMMS, 0, 0, 'r', 'LineWidth', 6);
rectangle('Position', [[0, 0] - 1, 2*1, 2*1]*2/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'r','LineStyle','--',LineWidth=2);
ylabel(col,'$N_{count}$','interpreter','latex')
col.TickLabelInterpreter = "latex";

savefig_FC([folderout filesep 'Urel_count'],8,6,'pdf')
savefig_FC([folderout filesep 'Urel_count'],8,6,'fig')
