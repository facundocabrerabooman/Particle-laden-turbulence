clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'TrCer_1000_13_fullg_inertial';
fname_tracers = 'trajs_TrCer_1000_13_fullg_tracers';

folderin = '/Volumes/landau1/TrCer_1000/dat/';
folderin_tracers = '/Volumes/landau1/Tracers/dat/';

folderout = folderin;
cd(folderin)

Fs=2990; % Frame rate

%% Import data

d = dat_to_mat(folderin, fname);
%save(fname,'d','-v7.3')

%% Track Particles (i.e. go from particle positions to trajectories)

maxdist = 1;  
lmin=10;
flag_pred=0;
npriormax=4;
porder=3;
flag_conf=1;
numFrames = 9e9;

[traj,~]=track3d_fc_stb(d,folderin,folderout,fname,maxdist,lmin,flag_pred,npriormax,porder,flag_conf, numFrames, Fs);

%% Only keep long tracks -- redundant if using track3d_fc_stb.m
L = arrayfun(@(X)(numel(X.x)),traj);
Ilong = find(L>=10);
%% Find proper filter width
if pi==pi
[s(1), m(1), w]=findFilterWidth_PTV(traj(Ilong),'x');
[s(2), m(2), w]=findFilterWidth_PTV(traj(Ilong),'y');
[s(3), m(3), w]=findFilterWidth_PTV(traj(Ilong),'z');
%% Set cool colors for plots
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

figure;
yyaxis left
loglog(w,s(1).vx,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on;
loglog(w,s(2).vx,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(w,s(3).vx,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);
hold off

yyaxis right
loglog(w,s(1).ax,'^-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on;
loglog(w,s(2).ax,'^-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(w,s(3).ax,'^-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

plot([10 10],ylim,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
yyaxis left
legend('$V_x$','$V_y$','$V_z$','$A_x$','$A_y$','$A_z$','interpreter','latex',Location='southwest',FontSize=12);
title('$std.(w)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$fliter\ width\ w$','interpreter','latex',FontWeight='bold',FontSize=18)

yyaxis left
ylabel('$\sigma_{v}$','interpreter','latex',FontWeight='bold',FontSize=24)
yyaxis right
ylabel('$\sigma_{a}$','interpreter','latex',FontWeight='bold',FontSize=24)

grid on
axis padded

folderout = ['filter_check_' fname filesep];
mkdir(folderout)
savefig_custom([folderout 'filter_check_' fname],8,6,'pdf')
savefig_custom([folderout 'filter_check_' fname],8,6,'fig')
save(['filter_check_' fname filesep 'output_filtering.mat'],'s','m','w')
end
%%  Estimate filtered tracks, velocities and accelerations with optimal filter width
wopt = 10;
lopt = 30;

Fs = 2990;

%[~, tracklong]=compute_vel_acc_traj(traj(Ilong),Fs,wopt,lopt);

tracklong=calcVelLEM(traj,wopt,lopt,Fs, wopt, lopt);

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong)==1);

save(['trajs_' fname '.mat'],'Ine','tracklong')

%% Traj. graph

clear color c
xt=vertcat(tracklong.x);
yt=vertcat(tracklong.y);
zt=vertcat(tracklong.z);

figure(10);clf

plot3(xt,yt,zt,'o');hold on

axis equal
box on
xlabel('mm')
ylabel('mm')
zlabel('mm')

%% Compute slip velocity

tracklong_inertial = tracklong; clear tracklong

load([folderin_tracers filesep fname_tracers '.mat'])
tracklong_tracers = tracklong; clear tracklong



stop



% Initialize arrays to store results
bigParticleVelocityCorrected = zeros(size(bigParticleTrajectory));
averageDustParticleVelocity = zeros(size(bigParticleTrajectory, 2), 3);

% Main loop for each time step
for t = 1:size(bigParticleTrajectory, 2)
    % Extract big particle position for the current time step
    bigParticlePosition = [bigParticleData(1).x(t); bigParticleData(1).y(t); bigParticleData(1).z(t)];

    % Extract dust particle positions at the current time step
    dustParticlePositions = zeros(3, numParticles);
    for particleIdx = 1:numParticles
        dustParticlePositions(:, particleIdx) = [dustParticleData(particleIdx).x(t); dustParticleData(particleIdx).y(t); dustParticleData(particleIdx).z(t)];
    end

    % Compute distances between big particle and dust particles
    distances = sqrt(sum((dustParticlePositions - bigParticlePosition).^2, 1));

    % Identify dust particles within the specified sphere
    particlesInSphere = distances < sphereDiameter / 2;

    % Compute average velocity of dust particles within the sphere
    if any(particlesInSphere)
        averageDustVelocity = mean(dustParticlePositions(:, particlesInSphere), 2);
    else
        averageDustVelocity = zeros(3, 1); % Default to zero if no particles in the sphere
    end

    % Store the average velocity for later use
    averageDustParticleVelocity(t, :) = averageDustVelocity';

    % Subtract average dust particle velocity from big particle velocity
    bigParticleVelocityCorrected(:, t) = bigParticleTrajectory(:, t) - averageDustVelocity;
end

% Display or use the results as needed
disp('Average Dust Particle Velocity:');
disp(averageDustParticleVelocity);

% The corrected big particle trajectory is stored in bigParticleVelocityCorrected















