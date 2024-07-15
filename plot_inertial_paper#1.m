%%  MAKE SPHERE PLOT
clfc

addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/'));

folderout = '/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/';
mkdir(folderout)
cd(folderout)

%%%%%%%%%%%%%%%%%%%%%% Plot

% % Find longest trajectories and their index.
% lengths_fullg = [];
% for i = 1:numel(trajs_conc_with_mean_field_fullg)
%     if isfield(trajs_conc_with_mean_field_fullg(i), 'Xf')
%         lengths_fullg(i,1) = length(trajs_conc_with_mean_field_fullg(i).Xf);
%         lengths_fullg(i,2) = i;
%     end
% end
% lengths_fullg = sortrows(lengths_fullg, 1, 'descend');
% %
% lengths_ddt = [];
% for i = 1:numel(trajs_conc_with_mean_field_ddt)
%     if isfield(trajs_conc_with_mean_field_ddt(i), 'Xf')
%         lengths_ddt(i,1) = length(trajs_conc_with_mean_field_ddt(i).Xf);
%         lengths_ddt(i,2) = i;
%     end
% end
% lengths_ddt = sortrows(lengths_ddt, 1, 'descend');
% %

% %%
% clc
% figure(1); clf; hold on
% term_vel_no_turb = 360;
%
% step = 10;
% for oo=1:2
% %vel_fullg = trajs_conc_with_mean_field_fullg(lengths_fullg(oo,2)).Vz(1:step:end);
% %scatter3(trajs_conc_with_mean_field_fullg(lengths_fullg(oo,2)).Xf(1:step:end),trajs_conc_with_mean_field_fullg(lengths_fullg(oo,2)).Yf(1:step:end),trajs_conc_with_mean_field_fullg(lengths_fullg(oo,2)).Zf(1:step:end),2,vel_fullg,'filled','Marker','o')
%
% %vel_ddt = sqrt((trajs_conc_with_mean_field_ddt(lengths_ddt(oo,2)).Vx(1:step:end)).^2 + (trajs_conc_with_mean_field_ddt(lengths_ddt(oo,2)).Vy(1:step:end)).^2 + (trajs_conc_with_mean_field_ddt(lengths_ddt(oo,2)).Vz(1:step:end)).^2);
% vel_ddt = abs(trajs_conc_with_mean_field_ddt(lengths_ddt(oo,2)).Vz(1:step:end));
% vel_ddt = vel_ddt./term_vel_no_turb;
% scatter3(trajs_conc_with_mean_field_ddt(lengths_ddt(oo,2)).Xf(1:step:end),trajs_conc_with_mean_field_ddt(lengths_ddt(oo,2)).Yf(1:step:end),trajs_conc_with_mean_field_ddt(lengths_ddt(oo,2)).Zf(1:step:end),50,vel_ddt,'filled','Marker','o')%'MarkerEdgeColor','white','LineWidth', 0.1)
% end
%
% colormap("autumn"); colorbar
%
% box on
% grid on
% xlabel('x'); xlim([-20 20])
% ylabel('y')
% zlabel('z'); zlim([-25 10])
% axis equal
%
% %%
% clc
% figure(1); clf; hold on
% term_vel_no_turb = 360;
%
% step = 10;
% for oo=1:2e1
% %vel_fullg = trajs_conc_without_mean_field_fullg(lengths_fullg(oo,2)).Vz(1:step:end);
% %scatter3(trajs_conc_without_mean_field_fullg(lengths_fullg(oo,2)).Xf(1:step:end),trajs_conc_without_mean_field_fullg(lengths_fullg(oo,2)).Yf(1:step:end),trajs_conc_without_mean_field_fullg(lengths_fullg(oo,2)).Zf(1:step:end),2,vel_fullg,'filled','Marker','o')
%
% %vel_ddt = sqrt((trajs_conc_without_mean_field_ddt(lengths_ddt(oo,2)).Vx(1:step:end)).^2 + (trajs_conc_without_mean_field_ddt(lengths_ddt(oo,2)).Vy(1:step:end)).^2 + (trajs_conc_without_mean_field_ddt(lengths_ddt(oo,2)).Vz(1:step:end)).^2);
% vel_ddt = abs(trajs_conc_without_mean_field_ddt(lengths_ddt(oo,2)).Vz(1:step:end));
% vel_ddt = vel_ddt./term_vel_no_turb;
% scatter3(trajs_conc_without_mean_field_ddt(lengths_ddt(oo,2)).Xf(1:step:end),trajs_conc_without_mean_field_ddt(lengths_ddt(oo,2)).Yf(1:step:end),trajs_conc_without_mean_field_ddt(lengths_ddt(oo,2)).Zf(1:step:end),50,vel_ddt,'filled','Marker','o')%'MarkerEdgeColor','white','LineWidth', 0.1)
% end
%
% colormap("autumn"); colorbar
%
% box on
% grid on
% xlabel('x'); xlim([-20 20])
% ylabel('y')
% zlabel('z'); zlim([-25 10])
% axis equal
%


%%% Concatenate data

trajs_conc_inertial = [];
load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/particle/ddt/trajsf_TrCer_1000_13_ddt_particle.mat')
trajs_conc_inertial = [trajs_conc_inertial tracklong];
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_inertial)==1);
clear tracklong
trajs_conc_inertial = trajs_conc_inertial(Ine);

trajs_conc_tracers = [];
load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/tracers/ddt/trajsf_TrCer_1000_13_ddt_tracer.mat')
trajs_conc_tracers = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_tracers)==1);
trajs_conc_tracers = trajs_conc_tracers(Ine);

%%%
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

%%
%%%
% Assuming you have two structures: trajs_conc_inertial and trajs_conc_tracers
% Parameters for the sphere
sphereRadius = 5;
%%%
for i=1:numel(trajs_conc_inertial)
    a(i,1) = mean(trajs_conc_inertial(i).Xf);
    a(i,2) = mean(trajs_conc_inertial(i).Yf);
    a(i,3) = mean(trajs_conc_inertial(i).Zf);
end
%%%
%%% MAKE SPHERE PLOT

inertialIdx = 192%:length(trajs_conc_inertial)
% Extract inertial particle information
Xf = trajs_conc_inertial(inertialIdx).Xf;
Zf = trajs_conc_inertial(inertialIdx).Yf;
Yf = trajs_conc_inertial(inertialIdx).Zf;
Tf = trajs_conc_inertial(inertialIdx).Tf;

timeIdx = numel(Tf)%round(numel(Tf)/2); % loop over time
time_inertial = Tf(timeIdx);

% Plot inertial particle
% Create a figure for the 3D plot
f = figure(10);clf
%f.Visible = 'off';
%f=figure;hold on
plot3(Xf(timeIdx), Yf(timeIdx), Zf(timeIdx), 'o', 'MarkerSize', 20,'MarkerFaceColor',color1);
box on
hold on;

% Loop over each element in trajs_conc_tracers
for tracerIdx = 1:length(trajs_conc_tracers)
    disp([sprintf('%.1f', inertialIdx/length(trajs_conc_inertial)*100) '--' sprintf('%.1f',tracerIdx/length(trajs_conc_tracers))])
    % Extract tracer particle information
    Xf_tracer = trajs_conc_tracers(tracerIdx).Xf;
    Zf_tracer = trajs_conc_tracers(tracerIdx).Yf;
    Yf_tracer = trajs_conc_tracers(tracerIdx).Zf;
    Tf_tracer = trajs_conc_tracers(tracerIdx).Tf;

    time_match=find(time_inertial == Tf_tracer);
    if ~isempty(time_match)
        % Check if the tracer particle is within the sphere for each point
        distance = sqrt((Xf(timeIdx) - Xf_tracer(time_match)).^2 + (Yf(timeIdx) - Yf_tracer(time_match)).^2 + (Zf(timeIdx) - Zf_tracer(time_match)).^2);
        flag_inside = find(distance <= sphereRadius);
        % Plot tracer particle if within the sphere
        if ~isempty(flag_inside) && (Zf_tracer(time_match) - Zf(timeIdx))>0
            % Loop over different view angles to create rotation effect
            plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'o', 'MarkerSize', 5,'MarkerFaceColor',mycolormap((size(mycolormap,1)+1)/2,:));hold on
        else
            plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'ko', 'MarkerSize', 2,'MarkerFaceColor','k');hold on
        end
    end
end

quiver3(Xf(timeIdx), Yf(timeIdx), Zf(timeIdx)+0.25, 0, 0, 3, 'LineWidth', 2, 'MaxHeadSize', 0.5,'LineWidth',4,'Color','red');

% % Plot the transparent sphere
[x, y, z] = sphere;
% h = surf( sphereRadius * x + Xf(timeIdx), sphereRadius * y + Yf(timeIdx), sphereRadius * z + Zf(timeIdx));
% alpha(h, 0.1);  % Set transparency (0 is completely transparent, 1 is opaque)

% Select only the upper half of the sphere (z >= 0)
x_half = x(z >= 0);
y_half = y(z >= 0);
z_half = z(z >= 0);

% Reshape the data to maintain the surface structure
x_half = reshape(x_half, [], size(x, 2));
y_half = reshape(y_half, [], size(y, 2));
z_half = reshape(z_half, [], size(z, 2));

h = surf(sphereRadius * x_half + Xf(timeIdx), sphereRadius * y_half + Yf(timeIdx), sphereRadius * z_half + Zf(timeIdx));
alpha(h, 0.15);  % Set transparency (0 is completely transparent, 1 is opaque)


view(-45, 10); % Set the view angle
% Set axis limits and labels
axis equal;
xlim([2,14]);
zlim([-10,0]);
ylim([-6,6]);

%xlabel('X (mm)');
%ylabel('Y (mm)');
%zlabel('Z (mm)');
%title(['Time: ' num2str(time_inertial)]);
xticks([])
yticks([])
zticks([])
savefig_FC('sphere',8,6,'fig')
savefig_FC('sphere',8,6,'pdf')

%%% MAKE VIDEO
if 1==pi
    % Set up VideoWriter
    videoFile = 'trajectory_video.mp4';  % Specify the desired file name
    writerObj = VideoWriter(videoFile, 'MPEG-4');
    writerObj.FrameRate = 5;  % Set the frame rate (adjust as needed)
    open(writerObj);

    % Loop over each element in trajs_conc_inertial
    for inertialIdx = 1:length(trajs_conc_inertial)
        % Extract inertial particle information
        Xf = trajs_conc_inertial(inertialIdx).Xf;
        Zf = trajs_conc_inertial(inertialIdx).Yf;
        Yf = trajs_conc_inertial(inertialIdx).Zf;
        Tf = trajs_conc_inertial(inertialIdx).Tf;


        for timeIdx = 1:5:numel(Tf) % loop over time
            time_inertial = Tf(timeIdx);

            % Plot inertial particle
            % Create a figure for the 3D plot
            f = figure(10); hold on
            f.Visible = 'off';
            %f=figure;hold on
            plot3(Xf(timeIdx), Yf(timeIdx), Zf(timeIdx), 'o', 'MarkerSize', 10,'MarkerFaceColor',color1);
            box on
            hold on;

            % Loop over each element in trajs_conc_tracers
            for tracerIdx = 1:length(trajs_conc_tracers)
                disp([sprintf('%.1f', inertialIdx/length(trajs_conc_inertial)*100) '--' sprintf('%.1f',tracerIdx/length(trajs_conc_tracers))])
                % Extract tracer particle information
                Xf_tracer = trajs_conc_tracers(tracerIdx).Xf;
                Zf_tracer = trajs_conc_tracers(tracerIdx).Yf;
                Yf_tracer = trajs_conc_tracers(tracerIdx).Zf;
                Tf_tracer = trajs_conc_tracers(tracerIdx).Tf;

                time_match=find(time_inertial == Tf_tracer);
                if ~isempty(time_match)
                    % Check if the tracer particle is within the sphere for each point
                    distance = sqrt((Xf(timeIdx) - Xf_tracer(time_match)).^2 + (Yf(timeIdx) - Yf_tracer(time_match)).^2 + (Zf(timeIdx) - Zf_tracer(time_match)).^2);

                    flag_inside = find(distance <= sphereRadius);
                    % Plot tracer particle if within the sphere
                    if ~isempty(flag_inside)
                        % Loop over different view angles to create rotation effect
                        plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'o', 'MarkerSize', 5,'MarkerFaceColor',mycolormap((size(mycolormap,1)+1)/2,:));hold on
                    else
                        plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'ko', 'MarkerSize', 1,'MarkerFaceColor','k');hold on
                    end
                end
            end

            % Plot the transparent sphere
            % [x, y, z] = sphere;
            % h = surf( sphereRadius * x + Xf(timeIdx), sphereRadius * y + Yf(timeIdx), sphereRadius * z + Zf(timeIdx));
            % alpha(h, 0.1);  % Set transparency (0 is completely transparent, 1 is opaque)

            view(-45, 10); % Set the view angle
            % Set axis limits and labels
            axis equal;
            xlim([-20,20]);
            zlim([-45,35]);
            ylim([-20,20]);
            xlabel('X-axis');
            ylabel('Y-axis');
            zlabel('Z-axis');
            title(['Time: ' num2str(time_inertial)]);
            writeVideo(writerObj, getframe(gcf));

            % Close the figure
            pause(1)
            close(gcf);

        end
        pause(1)

    end

    % Close the video file
    close(writerObj);
end


























































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SlipVel Stuff R6

clfc

addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/'));

folderout = '/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/inertial_cer_R6/g&ug/';
mkdir(folderout)
cd(folderout)

%folder_ddt = '/Volumes/landau1/TrCer_analysis_paper#1/exports/ddt_pairs/';
%folder_fullg = '/Volumes/landau1/TrCer_analysis_paper#1/exports/fullg_pairs/';
%folder_dec = '';

%%% Load slip velocities for fullg-ddt-dec

%load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/slipVeloData/slipVeloData_R_10/slipVelCylind_ddt_CONC.mat','AverSlipVelCylind_conc');
load('/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/inertial_cer_R6/data/slipVelCylind_ddt_CONC.mat','AverSlipVelCylind_conc');
AverSlipVelCylind_conc_ddt_R6 = AverSlipVelCylind_conc; clear AverSlipVelCylind_conc

%load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/slipVeloData/slipVeloData_R_10/slipVelCylind_fullg_CONC.mat','AverSlipVelCylind_conc');
load('/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/inertial_cer_R6/data/slipVelCylind_fullg_CONC.mat','AverSlipVelCylind_conc');
AverSlipVelCylind_conc_fullg_R6 = AverSlipVelCylind_conc; clear AverSlipVelCylind_conc


%%% Load particle velocities

% no turbulence case
load('/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/inertial_cer_R6/data/trajs_TrCer_1000_noturb.mat')
tracklong_noturb = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong_noturb)==1);
tracklong_noturb = tracklong_noturb(Ine);

%%% change reference system
counter=0;
tracklong_noturb_tmp = tracklong_noturb;
for i=1:numel(tracklong_noturb)
    tracklong_noturb(i).X = tracklong_noturb_tmp(i).x;
    tracklong_noturb(i).Y = tracklong_noturb_tmp(i).z;
    tracklong_noturb(i).Z = -tracklong_noturb_tmp(i).y;

    tracklong_noturb(i).Vx = tracklong_noturb_tmp(i).Vx;
    tracklong_noturb(i).Vy = tracklong_noturb_tmp(i).Vz;
    tracklong_noturb(i).Vz = -tracklong_noturb_tmp(i).Vy;

    tracklong_noturb(i).Ax = tracklong_noturb_tmp(i).Ax;
    tracklong_noturb(i).Ay = tracklong_noturb_tmp(i).Az;
    tracklong_noturb(i).Az = -tracklong_noturb_tmp(i).Ay;
end
clear tracklong_noturb_tmp


% Fullg particles -- in original reference system with y // to gravity
%load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/particle/trajsconc_fullg_particle.mat')
load('/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/inertial_cer_R6/data/trajsconc_fullg_particle.mat')
trajs_conc_with_mean_field_fullg = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_with_mean_field_fullg)==1);
trajs_conc_with_mean_field_fullg = trajs_conc_with_mean_field_fullg(Ine);

% DDT particles
%load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/particle/trajsconc_ddt_particle.mat')
load('/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/inertial_cer_R6/data/trajsconc_ddt_particle.mat')
trajs_conc_with_mean_field_ddt = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_with_mean_field_ddt)==1);
trajs_conc_with_mean_field_ddt = trajs_conc_with_mean_field_ddt(Ine);

% Dec particles
% load('/Volumes/landau1/TrCer_1000/dec/tracklongP_conc.mat')
% tracklong_dec = tracklong; clear tracklong


%%% Change reference system for particles

% FULLG
counter=0;
trajs_conc = trajs_conc_with_mean_field_fullg;

for i=1:numel(trajs_conc_with_mean_field_fullg)
    trajs_conc_with_mean_field_fullg(i).X = trajs_conc(i).x;
    trajs_conc_with_mean_field_fullg(i).Y = trajs_conc(i).z;
    trajs_conc_with_mean_field_fullg(i).Z = -trajs_conc(i).y;

    trajs_conc_with_mean_field_fullg(i).Vx = trajs_conc(i).Vx;
    trajs_conc_with_mean_field_fullg(i).Vy = trajs_conc(i).Vz;
    trajs_conc_with_mean_field_fullg(i).Vz = -trajs_conc(i).Vy;

    trajs_conc_with_mean_field_fullg(i).Ax = trajs_conc(i).Ax;
    trajs_conc_with_mean_field_fullg(i).Ay = trajs_conc(i).Az;
    trajs_conc_with_mean_field_fullg(i).Az = -trajs_conc(i).Ay;

    %trajs_conc(i).Z = trajs_conc(i).z;
end

clear trajs_conc

disp('GRAVITY goes downwards in z now; the vectors are anti-parallel')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DDT
counter=0;
trajs_conc = trajs_conc_with_mean_field_ddt;

for i=1:numel(trajs_conc_with_mean_field_ddt)
    trajs_conc_with_mean_field_ddt(i).X = trajs_conc(i).x;
    trajs_conc_with_mean_field_ddt(i).Y = trajs_conc(i).z;
    trajs_conc_with_mean_field_ddt(i).Z = -trajs_conc(i).y;

    trajs_conc_with_mean_field_ddt(i).Vx = trajs_conc(i).Vx;
    trajs_conc_with_mean_field_ddt(i).Vy = trajs_conc(i).Vz;
    trajs_conc_with_mean_field_ddt(i).Vz = -trajs_conc(i).Vy;

    trajs_conc_with_mean_field_ddt(i).Ax = trajs_conc(i).Ax;
    trajs_conc_with_mean_field_ddt(i).Ay = trajs_conc(i).Az;
    trajs_conc_with_mean_field_ddt(i).Az = -trajs_conc(i).Ay;

    %trajs_conc(i).Z = trajs_conc(i).z;
end

clear trajs_conc

disp('GRAVITY goes downwards in z now; the vectors are anti-parallel')


%%% Substract Mean Flow

% FULLG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Mean Velocity
trajs_conc_with_mean_field = trajs_conc_with_mean_field_fullg;

x = vertcat(trajs_conc_with_mean_field.Xf);
y = vertcat(trajs_conc_with_mean_field.Yf);
z = vertcat(trajs_conc_with_mean_field.Zf);

vx = vertcat(trajs_conc_with_mean_field.Vx);
vy = vertcat(trajs_conc_with_mean_field.Vy);
vz = vertcat(trajs_conc_with_mean_field.Vz);

numBinsX = 20;
numBinsY = 20;
numBinsZ = 30;

% Define the bin edges
xEdges = linspace(min(x), max(x), numBinsX+1);
yEdges = linspace(min(y), max(y), numBinsY+1);
zEdges = linspace(min(z), max(z), numBinsZ+1);

% Initialize matrices to store mean velocities
meanVx = zeros(numBinsX, numBinsY, numBinsZ);
meanVy = zeros(numBinsX, numBinsY, numBinsZ);
meanVz = zeros(numBinsX, numBinsY, numBinsZ);
count = zeros(numBinsX, numBinsY, numBinsZ); % To keep track of number of points in each bin

% Loop over bins and compute mean velocities
for i = 1:numBinsX
    for j = 1:numBinsY
        for k = 1:numBinsZ
            inBin = x >= xEdges(i) & x < xEdges(i+1) & ...
                y >= yEdges(j) & y < yEdges(j+1) & ...
                z >= zEdges(k) & z < zEdges(k+1);
            if any(inBin)
                meanVx(i, j, k) = mean(vx(inBin));
                meanVy(i, j, k) = mean(vy(inBin));
                meanVz(i, j, k) = mean(vz(inBin));
                count(i, j, k) = sum(inBin); % For debugging purposes
            end
        end
    end
end

% Generate bin centers for the interpolants
xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
zCenters = (zEdges(1:end-1) + zEdges(2:end)) / 2;

% Create griddedInterpolant objects for mean velocities
[Xc, Yc, Zc] = ndgrid(xCenters, yCenters, zCenters);
Fx = griddedInterpolant(Xc, Yc, Zc, meanVx, 'linear', 'nearest');
Fy = griddedInterpolant(Xc, Yc, Zc, meanVy, 'linear', 'nearest');
Fz = griddedInterpolant(Xc, Yc, Zc, meanVz, 'linear', 'nearest');

% Low pass filter the interpolant

xGrid = linspace(min(x), max(x), 20);
yGrid = linspace(min(y), max(y), 20);
zGrid = linspace(min(z), max(z), 30);

[X, Y, Z] = ndgrid(xGrid, yGrid, zGrid);

VxGrid = Fx(X, Y, Z);
VyGrid = Fy(X, Y, Z);
VzGrid = Fz(X, Y, Z);

sigma = 1; % Standard deviation of the Gaussian filter
VxFiltered = smooth3(VxGrid, 'gaussian', [5 5 5], sigma); % [5 5 5] is the filter size
VyFiltered = smooth3(VyGrid, 'gaussian', [5 5 5], sigma);
VzFiltered = smooth3(VzGrid, 'gaussian', [5 5 5], sigma);

% VxFiltered = VxGrid;
% VyFiltered = VyGrid;
% VzFiltered = VzGrid;

clear Fx Fy Fz
Fx = scatteredInterpolant(X(:), Y(:), Z(:), VxFiltered(:), 'linear', 'none');
Fy = scatteredInterpolant(X(:), Y(:), Z(:), VyFiltered(:), 'linear', 'none');
Fz = scatteredInterpolant(X(:), Y(:), Z(:), VzFiltered(:), 'linear', 'none');

% Substact mean flow
trajs_conc_without_mean_field_fullg = trajs_conc_with_mean_field_fullg;
for j=1:numel(trajs_conc_with_mean_field_fullg)

    for k = 1:numel(trajs_conc_with_mean_field_fullg(j).Vx)
        Vx(k) = trajs_conc_with_mean_field_fullg(j).Vx(k) - Fx(trajs_conc_with_mean_field_fullg(j).Xf(k),trajs_conc_with_mean_field_fullg(j).Yf(k),trajs_conc_with_mean_field_fullg(j).Zf(k));
        Vy(k) = trajs_conc_with_mean_field_fullg(j).Vy(k) - Fy(trajs_conc_with_mean_field_fullg(j).Xf(k),trajs_conc_with_mean_field_fullg(j).Yf(k),trajs_conc_with_mean_field_fullg(j).Zf(k));
        Vz(k) = trajs_conc_with_mean_field_fullg(j).Vz(k) - Fz(trajs_conc_with_mean_field_fullg(j).Xf(k),trajs_conc_with_mean_field_fullg(j).Yf(k),trajs_conc_with_mean_field_fullg(j).Zf(k));
    end
    trajs_conc_without_mean_field_fullg(j).Vx = [];
    trajs_conc_without_mean_field_fullg(j).Vx = Vx';
    trajs_conc_without_mean_field_fullg(j).Vy = [];
    trajs_conc_without_mean_field_fullg(j).Vy = Vy';
    trajs_conc_without_mean_field_fullg(j).Vz = [];
    trajs_conc_without_mean_field_fullg(j).Vz = Vz';
    clear Vx Vy Vz

end

clear Fx Fy Fz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DDT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Mean Velocity
trajs_conc_with_mean_field = trajs_conc_with_mean_field_ddt;

x = vertcat(trajs_conc_with_mean_field.Xf);
y = vertcat(trajs_conc_with_mean_field.Yf);
z = vertcat(trajs_conc_with_mean_field.Zf);

vx = vertcat(trajs_conc_with_mean_field.Vx);
vy = vertcat(trajs_conc_with_mean_field.Vy);
vz = vertcat(trajs_conc_with_mean_field.Vz);

numBinsX = 20;
numBinsY = 20;
numBinsZ = 30;

% Define the bin edges
xEdges = linspace(min(x), max(x), numBinsX+1);
yEdges = linspace(min(y), max(y), numBinsY+1);
zEdges = linspace(min(z), max(z), numBinsZ+1);

% Initialize matrices to store mean velocities
meanVx = zeros(numBinsX, numBinsY, numBinsZ);
meanVy = zeros(numBinsX, numBinsY, numBinsZ);
meanVz = zeros(numBinsX, numBinsY, numBinsZ);
count = zeros(numBinsX, numBinsY, numBinsZ); % To keep track of number of points in each bin

% Loop over bins and compute mean velocities
for i = 1:numBinsX
    for j = 1:numBinsY
        for k = 1:numBinsZ
            inBin = x >= xEdges(i) & x < xEdges(i+1) & ...
                y >= yEdges(j) & y < yEdges(j+1) & ...
                z >= zEdges(k) & z < zEdges(k+1);
            if any(inBin)
                meanVx(i, j, k) = mean(vx(inBin));
                meanVy(i, j, k) = mean(vy(inBin));
                meanVz(i, j, k) = mean(vz(inBin));
                count(i, j, k) = sum(inBin); % For debugging purposes
            end
        end
    end
end

% Generate bin centers for the interpolants
xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
zCenters = (zEdges(1:end-1) + zEdges(2:end)) / 2;

% Create griddedInterpolant objects for mean velocities
[Xc, Yc, Zc] = ndgrid(xCenters, yCenters, zCenters);
Fx = griddedInterpolant(Xc, Yc, Zc, meanVx, 'linear', 'nearest');
Fy = griddedInterpolant(Xc, Yc, Zc, meanVy, 'linear', 'nearest');
Fz = griddedInterpolant(Xc, Yc, Zc, meanVz, 'linear', 'nearest');

% Low pass filter the interpolant

xGrid = linspace(min(x), max(x), 20);
yGrid = linspace(min(y), max(y), 20);
zGrid = linspace(min(z), max(z), 30);

[X, Y, Z] = ndgrid(xGrid, yGrid, zGrid);

VxGrid = Fx(X, Y, Z);
VyGrid = Fy(X, Y, Z);
VzGrid = Fz(X, Y, Z);

sigma = 1; % Standard deviation of the Gaussian filter
VxFiltered = smooth3(VxGrid, 'gaussian', [5 5 5], sigma); % [5 5 5] is the filter size
VyFiltered = smooth3(VyGrid, 'gaussian', [5 5 5], sigma);
VzFiltered = smooth3(VzGrid, 'gaussian', [5 5 5], sigma);

% VxFiltered = VxGrid;
% VyFiltered = VyGrid;
% VzFiltered = VzGrid;

clear Fx Fy Fz
Fx = scatteredInterpolant(X(:), Y(:), Z(:), VxFiltered(:), 'linear', 'none');
Fy = scatteredInterpolant(X(:), Y(:), Z(:), VyFiltered(:), 'linear', 'none');
Fz = scatteredInterpolant(X(:), Y(:), Z(:), VzFiltered(:), 'linear', 'none');

% Substact mean flow
trajs_conc_without_mean_field_ddt = trajs_conc_with_mean_field_ddt;
for j=1:numel(trajs_conc_with_mean_field_ddt)

    for k = 1:numel(trajs_conc_with_mean_field_ddt(j).Vx)
        Vx(k) = trajs_conc_with_mean_field_ddt(j).Vx(k) - Fx(trajs_conc_with_mean_field_ddt(j).Xf(k),trajs_conc_with_mean_field_ddt(j).Yf(k),trajs_conc_with_mean_field_ddt(j).Zf(k));
        Vy(k) = trajs_conc_with_mean_field_ddt(j).Vy(k) - Fy(trajs_conc_with_mean_field_ddt(j).Xf(k),trajs_conc_with_mean_field_ddt(j).Yf(k),trajs_conc_with_mean_field_ddt(j).Zf(k));
        Vz(k) = trajs_conc_with_mean_field_ddt(j).Vz(k) - Fz(trajs_conc_with_mean_field_ddt(j).Xf(k),trajs_conc_with_mean_field_ddt(j).Yf(k),trajs_conc_with_mean_field_ddt(j).Zf(k));
    end
    trajs_conc_without_mean_field_ddt(j).Vx = [];
    trajs_conc_without_mean_field_ddt(j).Vx = Vx';
    trajs_conc_without_mean_field_ddt(j).Vy = [];
    trajs_conc_without_mean_field_ddt(j).Vy = Vy';
    trajs_conc_without_mean_field_ddt(j).Vz = [];
    trajs_conc_without_mean_field_ddt(j).Vz = Vz';
    clear Vx Vy Vz

end

clear Fx Fy Fz

%%% Compute distributions

% Particle Data
pdf_Vvert_without_meanflow_fullg = mkpdf5(trajs_conc_without_mean_field_fullg,'Vz',100,10);
pdf_Vvert_without_meanflow_ddt = mkpdf5(trajs_conc_without_mean_field_ddt,'Vz',100,10);

pdfVvert_with_meanflow_fullg = mkpdf5(trajs_conc_with_mean_field_fullg,'Vz',100,10);
pdfVvert_with_meanflow_ddt = mkpdf5(trajs_conc_with_mean_field_ddt,'Vz',100,10);

% Slip Vel.
%%% Create variable with velocity to use mkpdf5
for j=1:numel(AverSlipVelCylind_conc_fullg_R6)
    AverSlipVelCylind_conc_fullg_R6(j).VerticalVel = AverSlipVelCylind_conc_fullg_R6(j).Urel(:,3);
end
for j=1:numel(AverSlipVelCylind_conc_ddt_R6)
    AverSlipVelCylind_conc_ddt_R6(j).VerticalVel = AverSlipVelCylind_conc_ddt_R6(j).Urel(:,3); disp('watch out these 3s')
end

pdfSV_fullg = mkpdf5(AverSlipVelCylind_conc_fullg_R6,'VerticalVel',100,10); % in the direction of the particle motion rather than vertical
1
pdfSV_ddt = mkpdf5(AverSlipVelCylind_conc_ddt_R6,'VerticalVel',100,10);
2

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

%% Compute fluctuations

mean_raw_fullg = pdfVvert_with_meanflow_fullg.mean
std_raw_fullg = pdfVvert_with_meanflow_fullg.std
skew_raw_fullg = pdfVvert_with_meanflow_fullg.skewness
kurt_raw_fullg = pdfVvert_with_meanflow_fullg.flatness
%%
mean_SV_fullg = pdfSV_fullg.mean
std_SV_fullg = pdfSV_fullg.std
skew_SV_fullg = pdfSV_fullg.skewness
kurt_SV_fullg = pdfSV_fullg.flatness
%%
mean_raw_ddt = pdfVvert_with_meanflow_ddt.mean
std_raw_ddt = pdfVvert_with_meanflow_ddt.std
skew_raw_ddt = pdfVvert_with_meanflow_ddt.skewness
kurt_raw_ddt = pdfVvert_with_meanflow_ddt.flatness
%%
mean_SV_ddt = pdfSV_ddt.mean
std_SV_ddt = pdfSV_ddt.std
skew_SV_ddt = pdfSV_ddt.skewness
kurt_SV_ddt = pdfSV_ddt.flatness
%%
mean_womf_fullg = pdf_Vvert_without_meanflow_fullg.mean
std_womf_fullg = pdf_Vvert_without_meanflow_fullg.std
skew_womf_fullg = pdf_Vvert_without_meanflow_fullg.skewness
kurt_womf_fullg = pdf_Vvert_without_meanflow_fullg.flatness

%%
mean_womf_ddt = pdf_Vvert_without_meanflow_ddt.mean
std_womf_ddt = pdf_Vvert_without_meanflow_ddt.std
skew_womf_ddt = pdf_Vvert_without_meanflow_ddt.skewness
kurt_womf_ddt = pdf_Vvert_without_meanflow_ddt.flatness

%%
% vel3d=[];
%
% for i=1:numel(AverSlipVelCylind_conc_ddt_R10)
%     vel3d(i,:) = AverSlipVelCylind_conc_ddt_R10(i).Urelmean;
% end
%
% sigmax = var(vel3d(:,1))
% sigmay = var(vel3d(:,3))
% sigmaz = var(-vel3d(:,2))
%
% %%% fullg
%
% vel3d_fullg=[];
%
% for i=1:numel(AverSlipVelCylind_conc_fullg_R10)
%     vel3d_fullg(i,:) = AverSlipVelCylind_conc_fullg_R10(i).Urelmean;
% end
%
% sigmax = var(vel3d_fullg(:,1))
% sigmay = var(vel3d_fullg(:,3))
% sigmaz = var(-vel3d_fullg(:,2))



%% Plot PDFs not centered in zero

figure(1); hold on; clf

%%% plot slip vel
semilogy(pdfSV_fullg.xpdf ,pdfSV_fullg.pdf,'r-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfSV_fullg.mean,'r',LineWidth=2)

semilogy(pdfSV_ddt.xpdf ,pdfSV_ddt.pdf,'g-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfSV_ddt.mean,'g',LineWidth=2)

%%% plot particle with mean flow
semilogy(pdfVvert_with_meanflow_fullg.xpdf ,pdfVvert_with_meanflow_fullg.pdf,'m-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfVvert_with_meanflow_fullg(1).mean,'m',LineWidth=2)

semilogy(pdfVvert_with_meanflow_ddt.xpdf ,pdfVvert_with_meanflow_ddt.pdf,'c-o',MarkerSize=5,LineWidth=2);
xline(pdfVvert_with_meanflow_ddt(1).mean,'c',LineWidth=2)


%%% plot particle vel without mean flow
semilogy(pdf_Vvert_without_meanflow_fullg(1).xpdf,pdf_Vvert_without_meanflow_fullg(1).pdf,'b-o',MarkerSize=5,LineWidth=2);hold on
xline(pdf_Vvert_without_meanflow_fullg(1).mean,'b',LineWidth=3)

semilogy(pdf_Vvert_without_meanflow_ddt.xpdf,pdf_Vvert_without_meanflow_ddt(1).pdf,'b-o',MarkerSize=5,LineWidth=2);
xline(pdf_Vvert_without_meanflow_ddt(1).mean,'b',LineWidth=3)

%xline(pdfVvert_with_meanflow_ddt(3).mean,'y',LineWidth=3)

%%% plot no turb vel
xline(mean(vertcat(tracklong_noturb(:).Vz)),'k',LineWidth=3)

grid on
box on
axis padded

legend('$SlipV-FULLG$','Mean','$SlipV-DDT$','Mean','FULLG','Mean','DDT','Mean','No Turbulence - full g','interpreter','latex',Location='northeast');
%ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
%xlabel('Vertical Velocity (mm/s)',FontWeight='bold')
%xticks(-800:100:800)



folderout_tmp = 'pdfs';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'PDF_v_notcentered0'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'PDF_v_notcentered0'],8,6,'fig')


%% Plot PDFs centered in zero

figure(2); hold on; clf

%%% plot slip vel
f1=semilogy(pdfSV_fullg.xpdfn,pdfSV_fullg.pdfn,'-*',LineWidth=2,MarkerSize=8,Color=color3(1,:));hold on
f2=semilogy(pdfSV_ddt.xpdfn ,pdfSV_ddt.pdfn,'-^',LineWidth=2,MarkerSize=8,Color=color3(2,:));hold on

%%% plot particle with mean flow
f3=semilogy(pdfVvert_with_meanflow_fullg.xpdfn, pdfVvert_with_meanflow_fullg.pdfn,'-pentagram',LineWidth=2,MarkerSize=8,Color=color3(1,:));hold on
f4=semilogy(pdfVvert_with_meanflow_ddt.xpdfn,pdfVvert_with_meanflow_ddt.pdfn,'-v',LineWidth=2,MarkerSize=8,Color=color3(2,:));hold on

%%% W/o mean flow
f5=semilogy(pdf_Vvert_without_meanflow_fullg.xpdfn, pdf_Vvert_without_meanflow_fullg.pdfn,'-x',LineWidth=2,MarkerSize=8,Color=color3(1,:));hold on
f6=semilogy(pdf_Vvert_without_meanflow_ddt.xpdfn,pdf_Vvert_without_meanflow_ddt.pdfn,'-d',LineWidth=2,MarkerSize=8,Color=color3(2,:));hold on

xpdf=linspace(-5,5,1024);
f7=plot(xpdf,normpdf(xpdf,0,1),'k',LineWidth=3);

%%% plot particle vel without mean flow
%semilogy(pdf_Vvert_without_meanflow_fullg.xpdfn,pdf_Vvert_without_meanflow_fullg.pdfn,'b-o',MarkerSize=5,LineWidth=2);hold on
%semilogy(pdf_Vvert_without_meanflow_ddt.xpdfn,pdf_Vvert_without_meanflow_ddt.pdfn,'b-o',MarkerSize=5,LineWidth=2);

%%% plot no turb vel
%xline(mean(vertcat(tracklong_noturb(:).Vz)),'k',LineWidth=3)

grid on
box on
axis padded

legend('SlipV Earth','SlipV Microg.','Raw Earth','Raw  Microg.','W/o Mean Field Earth','W/o Mean Field Microg.','interpreter','latex',FontSize=20, Location = 'south')%, Position=[0.75, 0.75, 0.1, 0.2]); % Adjust as needed [left, bottom, width, height]);

ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
%xlabel('Vertical Velocity (mm/s)',FontWeight='bold')
%xticks(-800:100:800)
xlim([-7 7])

folderout_tmp = 'pdfs';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'PDF_v_centered0'],8,6,'pdf')

%%% INSET
%
% insetAxes = axes('Position', [0.5, 0.18, 0.05, 0.4]); % [left, bottom, width, height]
%
% %%% fullg-ddt-dec slip settling velocity
% mean2 = pdfSV_fullg.mean;
% std = pdfSV_fullg.std/2;
% h6 =  scatter(0.5,mean2,100,'MarkerFaceColor',color3(1,:),'MarkerEdgeColor',color3(1,:),'Marker','*'); hold on
%
% mean2 = pdfSV_ddt.mean;
% std = pdfSV_ddt.std/2;
% h7 =scatter(0.5,mean2,100,'MarkerFaceColor',color3(2,:),'MarkerEdgeColor',color3(2,:),'Marker','^'); hold on
%
% %%% fullg-ddt-dec particle settling velocity
% mean2 = pdfVvert_with_meanflow_fullg.mean;
% h1=scatter(0.5,mean2,100,'MarkerFaceColor',color3(1,:),'MarkerEdgeColor',color3(1,:),'Marker','pentagram'); hold on
%
% mean2 = pdfVvert_with_meanflow_ddt.mean;
% h2=scatter(0.5,mean2,100,'MarkerFaceColor',color3(2,:),'MarkerEdgeColor',color3(2,:),'Marker','v'); hold on
%
% mean2 = pdf_Vvert_without_meanflow_fullg.mean;
% h3=scatter(0.5,mean2,100,'MarkerFaceColor',color3(1,:),'MarkerEdgeColor',color3(1,:),'Marker','x'); hold on
%
% mean2= pdf_Vvert_without_meanflow_ddt.mean;
% h4=scatter(0.5,mean2,100,'MarkerFaceColor',color3(2,:),'MarkerEdgeColor',color3(2,:),'Marker','d'); hold on
%
% %%% plot no turb vel
% h5=scatter(0.5,mean(vertcat(tracklong_noturb.Vz)),100,'MarkerFaceColor','r');
%
%
% ylim([-400 400])
% xlim([0 1]); xticks([])
%legend([h1 h2 h3 h4 h5 h6 h7],'FULLG','DDT','FULLG without MEANFIELD','DDT without MEANFIELD','No Turbulence - full g','Slip V FULLG','Slip V DDT','interpreter','latex',Location='northeast');
%legend(h5,'Terminal Velocity in Earth','interpreter','latex','Location','southeast');%,Orientation='vertical');
%legend('show', 'Orientation', 'vertical', 'Location', 'northeast','Position', [0.95, 0.4, 0.1, 0.2]); % [left, bottom, width, height]

%ylabel('$\langle V \rangle$ (mm/s)','Interpreter','latex')
%box on; grid on



folderout_tmp = 'pdfs';
mkdir(folderout_tmp)
savefig_FC([folderout_tmp filesep 'PDF_v_centered0'],8,6,'pdf')
savefig_FC([folderout_tmp filesep 'PDF_v_centered0'],8,6,'fig')

%% Plot Slip Velocity Signal

figure(10);hold on; box on
vel3d=[];
start=1;
finish=100;

counter=0;
for i=start:finish
    counter=counter+1;
    vel3d(counter,:) = AverSlipVelCylind_conc_ddt_R6(i).Urelmean;
end

plot(vertcat(AverSlipVelCylind_conc_ddt_R6(start:finish).t)./2996,vel3d(:,1),'r.-')
%plot(vertcat(AverSlipVelCylind_conc_ddt_R10(start:finish).t),vel3d(:,1),'r.-')
plot(vertcat(AverSlipVelCylind_conc_ddt_R6(start:finish).t)./2996,vel3d(:,3),'g.-')
plot(vertcat(AverSlipVelCylind_conc_ddt_R6(start:finish).t)./2996,-vel3d(:,2),'b.-')

legend({'x','y','z (g)'})
xlabel('t (s)')
ylabel('Vel (mm/s)')
grid on;box on

savefig_FC('velocity_signal',8,6,'pdf')
savefig_FC('velocity_signal',8,6,'fig')



%% Plot Fluctuations and mean versus R

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

clc
figure(4);clf

Rtmp = [4 6 8 10 15:10:55];
for i=1:numel(Rtmp)
    i/numel(Rtmp)

    %try
    % fullg
    load(['/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/slipVeloData/slipVeloData_R_' num2str(Rtmp(i)) '/slipVelCylind_fullg_CONC.mat'],'AverSlipVelCylind_conc');
    for j=1:numel(AverSlipVelCylind_conc) % add new field to use mkfpdf5
        AverSlipVelCylind_conc(j).VerticalVel = AverSlipVelCylind_conc(j).Urel(:,3);
    end

    pdfSV_tmp_fullg = mkpdf5(AverSlipVelCylind_conc,'VerticalVel',100,10);


    clearvars -except i Rtmp color3 color1 pdfSV_tmp_fullg
    % ddt
    load(['/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/slipVeloData/slipVeloData_R_' num2str(Rtmp(i)) '/slipVelCylind_ddt_CONC.mat'],'AverSlipVelCylind_conc');
    for j=1:numel(AverSlipVelCylind_conc) % add new field to use mkfpdf5
        AverSlipVelCylind_conc(j).VerticalVel = AverSlipVelCylind_conc(j).Urel(:,3);
    end

    pdfSV_tmp_ddt = mkpdf5(AverSlipVelCylind_conc,'VerticalVel',100,10);

    %catch end


    %%% plot
    %%% mean

    figure(4);hold on
    subplot(2,1,1)
    scatter(Rtmp(i), pdfSV_tmp_fullg.mean, 'p', 'MarkerFaceColor', color3(1,:), 'MarkerEdgeColor', color3(1,:), 'SizeData', 300);
    hold on;
    scatter(Rtmp(i),pdfSV_tmp_ddt.mean,'^','MarkerFaceColor',color3(2,:),'MarkerEdgeColor',color3(2,:),'SizeData',300); hold on
    ylabel('$\langle \mathrm{Slip Vel.} \rangle$ (mm/s)','Interpreter','latex')
    xlim([0 60])
    ylim([150 350])
    xline(6,'r--','Linewidth',3)
    box on
    grid on


    %%% std
    subplot(2,1,2)
    scatter(Rtmp(i), pdfSV_tmp_fullg.std, 'p', 'MarkerFaceColor', color3(1,:), 'MarkerEdgeColor', color3(1,:), 'SizeData', 300);
    hold on;
    scatter(Rtmp(i),pdfSV_tmp_ddt.std,'^','MarkerFaceColor',color3(2,:),'MarkerEdgeColor',color3(2,:),'SizeData',300); hold on
    ylabel('$\mathrm{STD}$ (mm/s)','Interpreter','latex')
    xline(6,'r--','Linewidth',3)

end

box on
grid on
xlabel('R (mm)')
xlim([0 60])
ylim([165 200])
legend({'Earth','Microgravity'},FontSize=20,Location = 'southeast')

folderout_tmp = 'Different_Rs';
mkdir(folderout_tmp)
subplot(2,1,1)
savefig_FC([folderout_tmp filesep 'slip_versus_R'],8,6,'pdf')
%subplot(2,1,2)
savefig_FC([folderout_tmp filesep 'slip_versus_R'],8,6,'fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot 3D traj - graphic abstract

load('/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/inertial_cer/inertial_cer_R6/data/trajsconc_fullg_particle.mat')
%%%
figure(5)
ii_long = 629;
plot3(tracklong(ii_long).Xf,tracklong(ii_long).Yf,tracklong(ii_long).Zf,'.')

stop
ii_long=1;
for  ii=1:2052

    if ii_long<numel(tracklong(ii).Xf)
        ii_long = ii;
    end
end

ii_long

%% load trajs
trajs_conc_inertial = [];
load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/particle/fullg/trajsf_TrCer_1000_11_fullg_particle.mat')
trajs_conc_inertial = [trajs_conc_inertial tracklong];
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_inertial)==1);
clear tracklong
trajs_conc_inertial = trajs_conc_inertial(Ine);

% find longest trajectory to use as inertialIdx below
% ii_long=1;
% for  ii=1:numel(trajs_conc_inertial)
% 
% if ii_long<numel(trajs_conc_inertial(ii).Xf)
%     ii_long = ii;
% end
% end
% stop

trajs_conc_tracers = [];
load('/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports/tracers/fullg/trajsf_TrCer_1000_11_fullg_tracer.mat')
trajs_conc_tracers = tracklong; clear tracklong
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_tracers)==1);
trajs_conc_tracers = trajs_conc_tracers(Ine);

addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

%%
% Assuming you have two structures: trajs_conc_inertial and trajs_conc_tracers
% Parameters for the sphere

sphereRadius = 10;

%%%
for i=1:numel(trajs_conc_inertial)
    a(i,1) = mean(trajs_conc_inertial(i).Xf);
    a(i,2) = mean(trajs_conc_inertial(i).Yf);
    a(i,3) = mean(trajs_conc_inertial(i).Zf);
end
%%%
%%% MAKE SPHERE PLOT

inertialIdx = 68%:length(trajs_conc_inertial)
% Extract inertial particle information
Xf = trajs_conc_inertial(inertialIdx).Xf;
Zf = trajs_conc_inertial(inertialIdx).Yf;
Yf = trajs_conc_inertial(inertialIdx).Zf;
Tf = trajs_conc_inertial(inertialIdx).Tf;

f = figure(10);clf
counter=0;
for timeIdx = [1:20:80]%60:numel(Tf)%round(numel(Tf)/2); % loop over time
    time_inertial = Tf(timeIdx);
    counter = counter+1;
colorstr = ['r','g','b','k','m','c'];
    % Plot inertial particle
    % Create a figure for the 3D plot

    %f.Visible = 'off';
    %f=figure;hold on
    plot3(Xf(timeIdx), Yf(timeIdx), Zf(timeIdx), 'o', 'MarkerSize', 20,'MarkerFaceColor',color1); hold on
    box on
    hold on;

    % Loop over each element in trajs_conc_tracers
    for tracerIdx = 1:length(trajs_conc_tracers)
        disp([sprintf('%.1f', inertialIdx/length(trajs_conc_inertial)*100) '--' sprintf('%.1f',tracerIdx/length(trajs_conc_tracers))])
        % Extract tracer particle information
        Xf_tracer = trajs_conc_tracers(tracerIdx).Xf;
        Zf_tracer = trajs_conc_tracers(tracerIdx).Yf;
        Yf_tracer = trajs_conc_tracers(tracerIdx).Zf;
        Tf_tracer = trajs_conc_tracers(tracerIdx).Tf;
%        plot3(Xf_tracer, Yf_tracer, Zf_tracer,'k.')% [colorstr(counter) '-o'], 'MarkerFaceColor',colorstr(counter), 'MarkerSize', 10,'Filled');hold on


        time_match=find(time_inertial == Tf_tracer);
        if ~isempty(time_match)
            plot3(Xf_tracer, Yf_tracer, Zf_tracer,'k-.')% [colorstr(counter) '-o'], 'MarkerFaceColor',colorstr(counter), 'MarkerSize', 10,'Filled');hold on

            % Check if the tracer particle is within the sphere for each point
            distance = sqrt((Xf(timeIdx) - Xf_tracer(time_match)).^2 + (Yf(timeIdx) - Yf_tracer(time_match)).^2 + (Zf(timeIdx) - Zf_tracer(time_match)).^2);
            flag_inside = find(distance <= sphereRadius);
            % Plot tracer particle if within the sphere
            %if ~isempty(flag_inside) && (Zf_tracer(time_match) - Zf(timeIdx))>0
            if ~isempty(flag_inside) 
                plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'r-o', 'MarkerFaceColor','r', 'MarkerSize', 5);hold on
            %else
            %    plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'k.', 'MarkerSize', 1,'MarkerFaceColor','k');hold on
            end
        end
    end

    %quiver3(Xf(timeIdx), Yf(timeIdx), Zf(timeIdx)+0.25, 0, 0, 3, 'LineWidth', 2, 'MaxHeadSize', 0.5,'LineWidth',4,'Color','red'); hold on
end
xlim([-10 20])
ylim([-20 5])
zlim([-50 -20 ])
stop
% % Plot the transparent sphere
[x, y, z] = sphere;
% h = surf( sphereRadius * x + Xf(timeIdx), sphereRadius * y + Yf(timeIdx), sphereRadius * z + Zf(timeIdx));
% alpha(h, 0.1);  % Set transparency (0 is completely transparent, 1 is opaque)

% Select only the upper half of the sphere (z >= 0)
x_half = x(z >= 0);
y_half = y(z >= 0);
z_half = z(z >= 0);

% Reshape the data to maintain the surface structure
x_half = reshape(x_half, [], size(x, 2));
y_half = reshape(y_half, [], size(y, 2));
z_half = reshape(z_half, [], size(z, 2));

h = surf(sphereRadius * x_half + Xf(timeIdx), sphereRadius * y_half + Yf(timeIdx), sphereRadius * z_half + Zf(timeIdx));
alpha(h, 0.15);  % Set transparency (0 is completely transparent, 1 is opaque)


view(-45, 10); % Set the view angle
% Set axis limits and labels
axis equal;
xlim([2,14]);
zlim([-10,0]);
ylim([-6,6]);

%xlabel('X (mm)');
%ylabel('Y (mm)');
%zlabel('Z (mm)');
%title(['Time: ' num2str(time_inertial)]);
xticks([])
yticks([])
zticks([])
%savefig_FC('sphere',8,6,'fig')
%savefig_FC('sphere',8,6,'pdf')
stop

%%% MAKE VIDEO
if 1==pi
    % Set up VideoWriter
    videoFile = 'trajectory_video.mp4';  % Specify the desired file name
    writerObj = VideoWriter(videoFile, 'MPEG-4');
    writerObj.FrameRate = 5;  % Set the frame rate (adjust as needed)
    open(writerObj);

    % Loop over each element in trajs_conc_inertial
    for inertialIdx = 1:length(trajs_conc_inertial)
        % Extract inertial particle information
        Xf = trajs_conc_inertial(inertialIdx).Xf;
        Zf = trajs_conc_inertial(inertialIdx).Yf;
        Yf = trajs_conc_inertial(inertialIdx).Zf;
        Tf = trajs_conc_inertial(inertialIdx).Tf;


        for timeIdx = 1:5:numel(Tf) % loop over time
            time_inertial = Tf(timeIdx);

            % Plot inertial particle
            % Create a figure for the 3D plot
            f = figure(10); hold on
            f.Visible = 'off';
            %f=figure;hold on
            plot3(Xf(timeIdx), Yf(timeIdx), Zf(timeIdx), 'o', 'MarkerSize', 10,'MarkerFaceColor',color1);
            box on
            hold on;

            % Loop over each element in trajs_conc_tracers
            for tracerIdx = 1:length(trajs_conc_tracers)
                disp([sprintf('%.1f', inertialIdx/length(trajs_conc_inertial)*100) '--' sprintf('%.1f',tracerIdx/length(trajs_conc_tracers))])
                % Extract tracer particle information
                Xf_tracer = trajs_conc_tracers(tracerIdx).Xf;
                Zf_tracer = trajs_conc_tracers(tracerIdx).Yf;
                Yf_tracer = trajs_conc_tracers(tracerIdx).Zf;
                Tf_tracer = trajs_conc_tracers(tracerIdx).Tf;

                time_match=find(time_inertial == Tf_tracer);
                if ~isempty(time_match)
                    % Check if the tracer particle is within the sphere for each point
                    distance = sqrt((Xf(timeIdx) - Xf_tracer(time_match)).^2 + (Yf(timeIdx) - Yf_tracer(time_match)).^2 + (Zf(timeIdx) - Zf_tracer(time_match)).^2);

                    flag_inside = find(distance <= sphereRadius);
                    % Plot tracer particle if within the sphere
                    if ~isempty(flag_inside)
                        % Loop over different view angles to create rotation effect
                        plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'o', 'MarkerSize', 5,'MarkerFaceColor',mycolormap((size(mycolormap,1)+1)/2,:));hold on
                    else
                        plot3(Xf_tracer(time_match), Yf_tracer(time_match), Zf_tracer(time_match), 'ko', 'MarkerSize', 1,'MarkerFaceColor','k');hold on
                    end
                end
            end

            % Plot the transparent sphere
            % [x, y, z] = sphere;
            % h = surf( sphereRadius * x + Xf(timeIdx), sphereRadius * y + Yf(timeIdx), sphereRadius * z + Zf(timeIdx));
            % alpha(h, 0.1);  % Set transparency (0 is completely transparent, 1 is opaque)

            view(-45, 10); % Set the view angle
            % Set axis limits and labels
            axis equal;
            xlim([-20,20]);
            zlim([-45,35]);
            ylim([-20,20]);
            xlabel('X-axis');
            ylabel('Y-axis');
            zlabel('Z-axis');
            title(['Time: ' num2str(time_inertial)]);
            writeVideo(writerObj, getframe(gcf));

            % Close the figure
            pause(1)
            close(gcf);

        end
        pause(1)

    end

    % Close the video file
    close(writerObj);
end








