clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

folderin_fullg = '/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/exports/fullg/';
folderin_ddt = '/Volumes/landau1/TrCer_analysis_paper#1/tracer_analysis/exports/ddt/';

%folderout=folderin;
folderout = '/Users/fcb/Library/CloudStorage/OneDrive-Stanford/Manuscripts/rsi_ddt/tracers';
mkdir(folderout)
cd(folderout)

Fs=2996; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Concatenate data fullg
trajs_conc = [];

load([folderin_fullg filesep 'trajsf_Tracers_01_fullg.mat'],'tracklong')
trajs_conc = [trajs_conc tracklong];
1
load([folderin_fullg filesep 'trajsf_Tracers_04_fullg.mat'],'tracklong')
trajs_conc = [trajs_conc tracklong];
2

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
trajs_conc_tmp = trajs_conc;
for i=1:numel(trajs_conc)

        trajs_conc(i).x = trajs_conc_tmp(i).x;
        trajs_conc(i).y = trajs_conc_tmp(i).z;
        trajs_conc(i).z = -trajs_conc_tmp(i).y;

        trajs_conc(i).Xf = trajs_conc_tmp(i).Xf;
        trajs_conc(i).Yf = trajs_conc_tmp(i).Zf;
        trajs_conc(i).Zf = -trajs_conc_tmp(i).Yf;

        trajs_conc(i).Vx = trajs_conc_tmp(i).Vx;
        trajs_conc(i).Vy = trajs_conc_tmp(i).Vz;
        trajs_conc(i).Vz = -trajs_conc_tmp(i).Vy;

        trajs_conc(i).Ax = trajs_conc_tmp(i).Ax;
        trajs_conc(i).Ay = trajs_conc_tmp(i).Az;
        trajs_conc(i).Az = -trajs_conc_tmp(i).Ay;

        %trajs_conc(i).Z = trajs_conc(i).z;
end
clear trajs_conc_tmp
disp('GRAVITY goes downwards in z now; the vectors are anti-parallel')


%trajs_conc_with_mean_field = trajs_conc;

% try
% save([folderin filesep 'traj_conc'],'trajs_conc','-v7.3')
% catch
% end

% Keep only trajectories in the center of FOV
if pi==1
disp('Keep only trajectories in the center of FOV')
radius = 10;
center = [0,0,5];

filtered_trajs = struct('ntraj', {}, 'length', {}, 'Ntrack', {}, 't', {}, ...
                        'x', {}, 'y', {}, 'z', {}, 'nmatch', {}, ...
                        'Xf', {}, 'Yf', {}, 'Zf', {}, 'Vx', {}, ...
                        'Vy', {}, 'Vz', {}, 'Ax', {}, 'Ay', {}, ...
                        'Az', {}, 'Tf', {}, 'Tf_acc', {}, 'Ntrackf', {}, ...
                       'Ntrack_acc', {}, 'Vabs', {}, 'Aabs', {});


for i = 1:length(trajs_conc)
    x = trajs_conc(i).x;
    y = trajs_conc(i).y;
    z = trajs_conc(i).z;
    
    distances = sqrt((x - center(1)).^2 + (y - center(2)).^2 + (z - center(3)).^2);
    
    if any(distances <= radius)
        filtered_trajs(end+1) = trajs_conc(i); 
    end
end
end
%filtered_trajs = trajs_conc;

trajs_conc_with_mean_field = filtered_trajs; clear trajs_conc 

clearvars -except folderin_fullg folderout fname Fs mycolormap color3 color1 trajs_conc_with_mean_field folderin_ddt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Mean Velocity

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
trajs_conc_without_mean_field = trajs_conc_with_mean_field;
for j=1:numel(trajs_conc_with_mean_field)

    for k = 1:numel(trajs_conc_with_mean_field(j).Vx)
        Vx(k) = trajs_conc_with_mean_field(j).Vx(k) - Fx(trajs_conc_with_mean_field(j).Xf(k),trajs_conc_with_mean_field(j).Yf(k),trajs_conc_with_mean_field(j).Zf(k)); 
        Vy(k) = trajs_conc_with_mean_field(j).Vy(k) - Fy(trajs_conc_with_mean_field(j).Xf(k),trajs_conc_with_mean_field(j).Yf(k),trajs_conc_with_mean_field(j).Zf(k)); 
        Vz(k) = trajs_conc_with_mean_field(j).Vz(k) - Fz(trajs_conc_with_mean_field(j).Xf(k),trajs_conc_with_mean_field(j).Yf(k),trajs_conc_with_mean_field(j).Zf(k)); 
    end
    trajs_conc_without_mean_field(j).Vx = [];
    trajs_conc_without_mean_field(j).Vx = Vx';
    trajs_conc_without_mean_field(j).Vy = [];
    trajs_conc_without_mean_field(j).Vy = Vy';
    trajs_conc_without_mean_field(j).Vz = [];
    trajs_conc_without_mean_field(j).Vz = Vz';
    clear Vx Vy Vz

end


trajs_conc_without_mean_field_fullg = trajs_conc_without_mean_field;
trajs_conc_with_mean_field_fullg = trajs_conc_with_mean_field;

clear trajs_conc_without_mean_field trajs_conc_with_mean_field

%% Concatenate data DDT
trajs_conc = [];

load([folderin_ddt filesep 'trajsf_Tracers_01_ddt.mat'],'tracklong')
trajs_conc = [trajs_conc tracklong];
1
load([folderin_ddt filesep 'trajsf_Tracers_02_ddt.mat'],'tracklong')
trajs_conc = [trajs_conc tracklong];

load([folderin_ddt filesep 'trajsf_Tracers_03_ddt.mat'],'tracklong')
trajs_conc = [trajs_conc tracklong];

load([folderin_ddt filesep 'trajsf_Tracers_04_ddt.mat'],'tracklong')
trajs_conc = [trajs_conc tracklong];
2

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
trajs_conc_tmp = trajs_conc;
counter=0;
for i=1:numel(trajs_conc)

        trajs_conc(i).x = trajs_conc_tmp(i).x;
        trajs_conc(i).y = trajs_conc_tmp(i).z;
        trajs_conc(i).z = -trajs_conc_tmp(i).y;

        trajs_conc(i).Xf = trajs_conc_tmp(i).Xf;
        trajs_conc(i).Yf = trajs_conc_tmp(i).Zf;
        trajs_conc(i).Zf = -trajs_conc_tmp(i).Yf;

        trajs_conc(i).Vx = trajs_conc_tmp(i).Vx;
        trajs_conc(i).Vy = trajs_conc_tmp(i).Vz;
        trajs_conc(i).Vz = -trajs_conc_tmp(i).Vy;

        trajs_conc(i).Ax = trajs_conc_tmp(i).Ax;
        trajs_conc(i).Ay = trajs_conc_tmp(i).Az;
        trajs_conc(i).Az = -trajs_conc_tmp(i).Ay;

        %trajs_conc(i).Z = trajs_conc(i).z;
end
clear trajs_conc_tmp


% Keep only trajectories in the center of FOV
if pi==1
disp('Keep only trajectories in the center of FOV')
radius = 10;
center = [0,0,5];

filtered_trajs = struct('ntraj', {}, 'length', {}, 'Ntrack', {}, 't', {}, ...
                        'x', {}, 'y', {}, 'z', {}, 'nmatch', {}, ...
                        'Xf', {}, 'Yf', {}, 'Zf', {}, 'Vx', {}, ...
                        'Vy', {}, 'Vz', {}, 'Ax', {}, 'Ay', {}, ...
                        'Az', {}, 'Tf', {}, 'Tf_acc', {}, 'Ntrackf', {}, ...
                       'Ntrack_acc', {}, 'Vabs', {}, 'Aabs', {});


for i = 1:length(trajs_conc)
    x = trajs_conc(i).x;
    y = trajs_conc(i).y;
    z = trajs_conc(i).z;
    
    distances = sqrt((x - center(1)).^2 + (y - center(2)).^2 + (z - center(3)).^2);
    
    if any(distances <= radius)
        filtered_trajs(end+1) = trajs_conc(i); 
    end
end
end
filtered_trajs = trajs_conc;

trajs_conc_with_mean_field = filtered_trajs; clear trajs_conc 

clearvars -except trajs_conc_with_mean_field_fullg trajs_conc_without_mean_field_fullg filtered_trajs folderin_ddt folderin_fullg folderout fname Fs mycolormap color3 color1 trajs_conc_with_mean_field 

%%% If doing microgravity data Only keep data before 5400frames, where the peak in acceleration
%%% happens +
%disp('microgravity?'); pause
% index = [];
% for i=1:numel(trajs_conc_with_mean_field)
%     if trajs_conc_with_mean_field(i).t(end)<0.5
%         index = [index i];
%     end
% end
% %trajs_conc_or = trajs_conc;
% trajs_conc_with_mean_field(index)=[];

disp('GRAVITY goes downwards in z now; the vectors are anti-parallel')


%clearvars -except folderin folderout fname Fs mycolormap color3 color1 trajs_conc trajs_conc_with_mean_field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Mean Velocity

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
Fx = griddedInterpolant(Xc, Yc, Zc, meanVx, 'cubic', 'nearest');
Fy = griddedInterpolant(Xc, Yc, Zc, meanVy, 'cubic', 'nearest');
Fz = griddedInterpolant(Xc, Yc, Zc, meanVz, 'cubic', 'nearest');

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

clear Fx Fy Fz
Fx = scatteredInterpolant(X(:), Y(:), Z(:), VxFiltered(:), 'linear', 'none');
Fy = scatteredInterpolant(X(:), Y(:), Z(:), VyFiltered(:), 'linear', 'none');
Fz = scatteredInterpolant(X(:), Y(:), Z(:), VzFiltered(:), 'linear', 'none');


% Substact mean flow
trajs_conc_without_mean_field = trajs_conc_with_mean_field;
for j=1:numel(trajs_conc_with_mean_field)

    for k = 1:numel(trajs_conc_with_mean_field(j).Vx)
        Vx(k) = trajs_conc_with_mean_field(j).Vx(k) - Fx(trajs_conc_with_mean_field(j).Xf(k),trajs_conc_with_mean_field(j).Yf(k),trajs_conc_with_mean_field(j).Zf(k)); 
        Vy(k) = trajs_conc_with_mean_field(j).Vy(k) - Fy(trajs_conc_with_mean_field(j).Xf(k),trajs_conc_with_mean_field(j).Yf(k),trajs_conc_with_mean_field(j).Zf(k)); 
        Vz(k) = trajs_conc_with_mean_field(j).Vz(k) - Fz(trajs_conc_with_mean_field(j).Xf(k),trajs_conc_with_mean_field(j).Yf(k),trajs_conc_with_mean_field(j).Zf(k)); 
    end
    trajs_conc_without_mean_field(j).Vx = [];
    trajs_conc_without_mean_field(j).Vx = Vx';
    trajs_conc_without_mean_field(j).Vy = [];
    trajs_conc_without_mean_field(j).Vy = Vy';
    trajs_conc_without_mean_field(j).Vz = [];
    trajs_conc_without_mean_field(j).Vz = Vz';
    clear Vx Vy Vz

end

trajs_conc_without_mean_field_ddt = trajs_conc_without_mean_field;
trajs_conc_with_mean_field_ddt = trajs_conc_with_mean_field;

clear trajs_conc_without_mean_field trajs_conc_with_mean_field trajs_conc trackout

%% Save trajs

%save('trajs','trajs_conc_with_mean_field_ddt','trajs_conc_without_mean_field_ddt','trajs_conc_with_mean_field_fullg','trajs_conc_without_mean_field_fullg','-v7.3')
%load('trajs','trajs_conc_with_mean_field_ddt','trajs_conc_without_mean_field_ddt','trajs_conc_with_mean_field_fullg','trajs_conc_without_mean_field_fullg')

%load('trajs','trajs_conc_with_mean_field_fullg')


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S2 - Epsilon - Eulerian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eulerian 2-point statistics -- without mean field
%% Compute S2 in intervals

for n = 0:10
    n
    step = 1e4;
start_idx = n * step + 1;
end_idx = min((n + 1) * step, length(trajs_conc_with_mean_field_fullg));
[eulerStats_tmp, pair_ddt] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_with_mean_field_fullg(start_idx:end_idx),[0.5 100],15,'off');
eulerStats_g{n+1} = eulerStats_tmp; clear eulerStats_ddt_tmp
end
eulerstats_g_fullg = eulerStats_g; clear eulerStats_g
save('eulerstats_fullg_withmf_fullgFOV','eulerstats_g_fullg');

for n = 0:10
    n
    step = 1e4;
start_idx = n * step + 1;
end_idx = min((n + 1) * step, length(trajs_conc_with_mean_field_ddt));
[eulerStats_tmp, pair_ddt] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_with_mean_field_ddt(start_idx:end_idx),[0.5 100],15,'off');
eulerStats_g{n+1} = eulerStats_tmp; clear eulerStats_ddt_tmp
end
eulerstats_g_ddt = eulerStats_g; clear eulerStats_g
save('eulerstats_ddt_withmf_fullgFOV','eulerstats_g_ddt');

%% Average subset S2s
load('eulerstats_fullg_withmf_fullgFOV.mat')

S2_aver=[];
for i=3:numel(eulerStats_g_fullg)
    S2_aver = [S2_aver eulerStats_g_fullg{i}.Splong{1,2}];
end

S2_aver_fullg = mean(S2_aver,2);
%%%
load('eulerstats_ddt_withmf_fullgFOV.mat')

S2_aver=[];
for i=3:numel(eulerStats_g_ddt)
    S2_aver = [S2_aver eulerStats_g_ddt{i}.Splong{1,2}];
end

S2_aver_ddt = mean(S2_aver,2);

% %% Plot S2 -- are they good? Testing
% for n=1:10
% figure
% 
% r = linspace(1,10,100);
% loglog(r,3e9/1e6*r.^(2/3),'--',Color=color3(3,:),LineWidth=2); hold on
% 
% loglog(eulerStats_g{n}.r,eulerStats_g{n}.Splong{1,2},'go-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
% 
% 
% end
% 
% stop

%% old
% Restricted to sphere 
%[eulerStats_fullg, pair_fullg] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_without_mean_field_fullg,[0.5 100],15,'off');
%[eulerStats_fullg, pair_fullg] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_without_mean_field_fullg,[0.5 100],15,'off');
%[eulerStats_ddt, pair_ddt] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_without_mean_field_ddt(1e4:end),[0.5 100],15,'off');

%% Plot S2 with epsilon inset
figure(3);clf

loglog(eulerstats_g_fullg{1}.r,S2_aver_fullg,'pentagram-',MarkerSize=8,Color=color3(1,:),LineWidth=2); hold on
loglog(eulerstats_g_ddt{1}.r,S2_aver_ddt,'^-',MarkerSize=8,Color=color3(2,:),LineWidth=2)
%loglog(eulerStats_ddt.r,eulerStats_ddt.Splong{1,2},'go-',MarkerSize=8,Color=color3(2,:),LineWidth=2);

r = linspace(1,10,100);
loglog(r,4e9/1e6*r.^(2/3),'--',Color=color3(3,:),LineWidth=2); hold on

ylabel('$S_2^{E,\parallel} (mm^2 s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded

legend('Earth','Microgravity','$r^{2/3}$','interpreter','latex',Location='northwest',FontSize=20)

%savefig_FC([folderout filesep 'S2_and_eps_womeanfield'],8,6,'pdf')% just to change font

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot epsilon
insetAxes = axes('Position', [0.55, 0.2, 0.3, 0.2]); % [left, bottom, width, height]

Ckolomogrov = 2.1;
%loglog(eulerStats_fullg.r,(eulerStats_fullg.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats_fullg.r','d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);
%hold on
%loglog(eulerStats_ddt.r,(eulerStats_ddt.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats_ddt.r','d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);

loglog(eulerstats_g_fullg{1}.r,(S2_aver_fullg./Ckolomogrov).^(1.5)./eulerstats_g_fullg{1}.r','pentagram-',MarkerSize=8,Color=color3(1,:),LineWidth=2); hold on
loglog(eulerstats_g_ddt{1}.r,(S2_aver_ddt./Ckolomogrov).^(1.5)./eulerstats_g_ddt{1}.r','^-',MarkerSize=8,Color=color3(2,:),LineWidth=2); hold on


r = linspace(1,10,100);
%%% compute epsilon
epsilon_fullg = (S2_aver_fullg./Ckolomogrov).^(1.5)./eulerstats_g_fullg{1}.r';
epsilon_fullg = mean(epsilon_fullg(4:7));

epsilon_ddt = (S2_aver_ddt./Ckolomogrov).^(1.5)./eulerstats_g_ddt{1}.r';
epsilon_ddt = mean(epsilon_ddt(4:7));
%%%
loglog(r,ones(size(r)).*epsilon_fullg,'--',Color=color3(1,:),LineWidth=2);hold on
loglog(r,ones(size(r)).*epsilon_ddt,'--',Color=color3(2,:),LineWidth=2)

set(gca,FontSize=20)
%ylabel('$\epsilon = (S_2^{E,\parallel}/C_k)^{3/2}\cdot r^{-1}$','interpreter','latex',FontWeight='bold',FontSize=20)
ylabel('$\epsilon (mm^2 s^{-3})$','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r (mm)$','interpreter','latex')
%ylim([9.5e3/1e6 3e4/1e6])
xlim([0.6 20])
xticks([1e0 1e1 2e1])
ylim([5e4 7e4])
legend('','','$\epsilon$','interpreter','latex',Location='northeast',FontSize=20)
%grid on
stop
savefig_FC([folderout filesep 'S2_and_eps_without_meanfield_fullfov_withmf'],8,6,'pdf')
savefig_FC([folderout filesep 'S2_and_eps_without_meanfield_fullfov_withmf'],8,6,'fig')

%% L int BATCHELOR PARAMETRIZATION GOOD!

load('eulerstats_fullg_withoutmeanfield_reducedFOV_r10.mat')
%%% Average subset S2s
S2_aver=[];
for i=3:numel(eulerStats_g_fullg)
    S2_aver = [S2_aver eulerStats_g_fullg{i}.Splong{1,2}];
end

S2_aver_fullg = mean(S2_aver,2);
%%%
load('eulerstats_ddt_withoutmeanfield_reducedFOV_r10.mat')

S2_aver=[];
for i=3:numel(eulerStats_g_ddt)
    S2_aver = [S2_aver eulerStats_g_ddt{i}.Splong{1,2}];
end

S2_aver_ddt = mean(S2_aver,2);
%%% 

figure(4);clf

%loglog(eulerStats_g_fullg{1}.r,S2_aver_fullg,'pentagram-',MarkerSize=8,Color=color3(1,:),LineWidth=2); hold on
loglog(eulerStats_g_ddt{1}.r(1:end-3),S2_aver_ddt(1:end-3),'^-',MarkerSize=8,Color=color3(2,:),LineWidth=2);hold on

ylabel('$S_2^{E,\parallel} (mm^2 s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded


%%%% FIT
r = eulerStats_g_ddt{1}.r'./1e3;
r=r(1:end-3);
%r = linspace(0,0.05,100)';

nu = 1e-6;
sigmaU = sqrt(1.6e4/1e6);
epsilon = 0.058;
zeta = 2/3;

[S2batch, eta, L, Lint, Ceps] = batchelorS2(r, nu, sigmaU, epsilon, zeta);

% Plot the results
loglog(r.*1e3, S2batch.*1e6, 'k-','LineWidth',2);

legend('Microgravity','Batchelor Parametrization','Interpreter','latex','Location','southeast')

savefig_FC([folderout filesep 'Lint'],8,6,'pdf')
savefig_FC([folderout filesep 'Lint'],8,6,'fig')
%% 
r = eulerStats_g_fullg{1}.r(1:end-3)';
y = S2_aver_fullg(1:end-3);
r=r./1e3;
y=y./1e6;

%% DOES NOT WORK Plot integral of R = 1-S2/sigma = Integral Scale - Without mean field
%load('trajs_with_without_mean_FOV_r10.mat')
load('eulerstats_ddt_withmf_fullgFOV.mat')
load('eulerstats_fullg_withmf_fullgFOV.mat')

figure(4);clf; hold on

Vx = vertcat(trajs_conc_with_mean_field_fullg.Vx);
Vy = vertcat(trajs_conc_with_mean_field_fullg.Vy);
Vz = vertcat(trajs_conc_with_mean_field_fullg.Vz);
sigma2 = (std(Vx).^2+std(Vy).^2+std(Vz).^2)/3;

%r = eulerStats_g_fullg{end}.Ruur;
%R = eulerStats_g_fullg{end}.Ruu;
%loglog(r,Splong,'k.-')
%stop
%R = 1-Splong/(27000);
%plot(r(1:end-3),R(1:end-3),'bo-') % plot correlation
%stop

r = eulerstats_g_fullg{1}.r';
Splong = S2_aver_fullg;
R = 1-Splong/(2*sigma2);
%semilogy(r,R,'.-')
%stop 
cumulative_integral = cumtrapz(r(1:end-3), R(1:end-3));
plot(r(1:end-3),cumulative_integral,'ro-')

%%%%%%%%%%%%%%% DDT Lint
%figure(5);clf
clear Vx Vy Vz r R Splong sigma2

Vx = vertcat(trajs_conc_with_mean_field_ddt.Vx);
Vy = vertcat(trajs_conc_with_mean_field_ddt.Vy);
Vz = vertcat(trajs_conc_with_mean_field_ddt.Vz);
sigma2 = (std(Vx).^2+std(Vy).^2+std(Vz).^2)/3;

%r = eulerStats_ddt.Ruur;
%R = eulerStats_ddt.Ruu;
%Splong = eulerStats_ddt.SplongAbs{1,2};
%loglog(r,Splong,'k.-')
%R = 1-eulerStats_ddt.SplongAbs{1,2}/(27000);
% plot(r,R,'bo-') % plot correlation
% stop

r = eulerstats_g_ddt{1}.r';
R = 1-S2_aver_ddt/(2*sigma2);
cumulative_integral = cumtrapz(r(1:end-3), R(1:end-3));
plot(r(1:end-3),cumulative_integral,'go-')

xlabel('r (mm)')
ylabel('Lint (mm)')

savefig_FC([folderout filesep 'int_scale' filesep 'L_int_fullFOV_withMF'],8,6,'pdf')
savefig_FC([folderout filesep 'int_scale' filesep 'L_int_fullFOV_withMF'],8,6,'fig')

%% Eulerian 2-point statistics -- WITH mean field
%% Compute S2 in intervals

for n = 0:10
    n
    step = 1e4;
start_idx = n * step + 1;
end_idx = min((n + 1) * step, length(trajs_conc_with_mean_field_fullg));
[eulerStats_tmp, pair_ddt] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_with_mean_field_fullg(start_idx:end_idx),[0.5 100],15,'off');
eulerStats_g{n+1} = eulerStats_tmp; clear eulerStats_ddt_tmp
end
eulerStats_g_fullg = eulerStats_g;
save('eulerstats_fullg_withmeanfield','eulerStats_g');
%%
for n = 0:10
    n
    step = 1e4;
start_idx = n * step + 1;
end_idx = min((n + 1) * step, length(trajs_conc_with_mean_field_ddt));
[eulerStats_tmp, pair_ddt] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_with_mean_field_ddt(start_idx:end_idx),[0.5 100],15,'off');
eulerStats_g{n+1} = eulerStats_tmp; clear eulerStats_ddt_tmp
end
eulerStats_g_ddt = eulerStats_g;
save('eulerstats_ddt_withmeanfield','eulerStats_g');

%% Average subset S2s
load('eulerstats_fullg_withmeanfield.mat')

S2_aver=[];
for i=3:numel(eulerStats_g_fullg)
    S2_aver = [S2_aver eulerStats_g_fullg{i}.Splong{1,2}];
end

S2_aver_fullg = mean(S2_aver,2);
%%%
load('eulerstats_ddt_withmeanfield.mat')

S2_aver=[];
for i=3:numel(eulerStats_g_ddt)
    S2_aver = [S2_aver eulerStats_g_ddt{i}.Splong{1,2}];
end

S2_aver_ddt = mean(S2_aver,2);

%% DOES NOT WORK Plot integral of R = 1-S2/sigma = Integral Scale - Without mean field
load('trajs_with_with_mean_reducedFOV_r10.mat')
load('eulerstats_ddt_withmeanfield.mat')
load('eulerstats_fullg_withmeanfield.mat')

figure(4);clf; hold on

Vx = vertcat(trajs_conc_with_mean_field_fullg.Vx);
Vy = vertcat(trajs_conc_with_mean_field_fullg.Vy);
Vz = vertcat(trajs_conc_with_mean_field_fullg.Vz);
sigma2 = (std(Vx).^2+std(Vy).^2+std(Vz).^2)/3;
stop
%r = eulerStats_g_fullg{end}.Ruur;
%R = eulerStats_g_fullg{end}.Ruu;
%loglog(r,Splong,'k.-')
%stop
%R = 1-Splong/(27000);
%plot(r(1:end-3),R(1:end-3),'bo-') % plot correlation
%stop

r = eulerStats_g_fullg{1}.r';
Splong = S2_aver_fullg;
%Splong = eulerStats_g_fullg{end}.Splong{1,2};
R = 1-Splong/(2*sigma2);

cumulative_integral = cumtrapz(r(1:end-3), R(1:end-3));
plot(r(1:end-3),cumulative_integral,'ro-')

%%%%%%%%%%%%%%% DDT Lint
%figure(5);clf
clear Vx Vy Vz r R Splong sigma2

Vx = vertcat(trajs_conc_without_mean_field_ddt.Vx);
Vy = vertcat(trajs_conc_without_mean_field_ddt.Vy);
Vz = vertcat(trajs_conc_without_mean_field_ddt.Vz);
sigma2 = (std(Vx).^2+std(Vy).^2+std(Vz).^2)/3;

%r = eulerStats_ddt.Ruur;
%R = eulerStats_ddt.Ruu;
%Splong = eulerStats_ddt.SplongAbs{1,2};
%loglog(r,Splong,'k.-')
%R = 1-eulerStats_ddt.SplongAbs{1,2}/(27000);
 %plot(r,R,'bo-') % plot correlation
 %stop

r = eulerStats_g_ddt{1}.r';
R = 1-S2_aver_ddt/(2*sigma2);
cumulative_integral = cumtrapz(r(1:end-3), R(1:end-3));

plot(r(1:end-3),cumulative_integral,'go-')

xlabel('r (mm)')
ylabel('Lint (mm)')


savefig_FC([folderout filesep 'int_scale' filesep 'L_int_withmf'],8,6,'pdf')
savefig_FC([folderout filesep 'int_scale' filesep 'L_int_withmf'],8,6,'fig')


%% Plot PDFs

% pdfV_fullg(1) = mkpdf5(trajs_conc_without_mean_field_fullg,'Vx',256,10);
% 1
% pdfV_fullg(2) = mkpdf5(trajs_conc_without_mean_field_fullg,'Vy',256,10);
% 2
% pdfV_fullg(3) = mkpdf5(trajs_conc_without_mean_field_fullg,'Vz',256,10);
% 3
% 
% pdfV_ddt(1) = mkpdf5(trajs_conc_without_mean_field_ddt,'Vx',256,10);
% 1
% pdfV_ddt(2) = mkpdf5(trajs_conc_without_mean_field_ddt,'Vy',256,10);
% 2
% pdfV_ddt(3) = mkpdf5(trajs_conc_without_mean_field_ddt,'Vz',256,10);
% 3
% 
% pdfA_fullg(1) = mkpdf5(trajs_conc_without_mean_field_fullg,'Ax',256,20);
% 4
% pdfA_fullg(2) = mkpdf5(trajs_conc_without_mean_field_fullg,'Ay',256,20);
% 5
% pdfA_fullg(3) = mkpdf5(trajs_conc_without_mean_field_fullg,'Az',256,20);
% 6
% 
% pdfA_ddt(1) = mkpdf5(trajs_conc_without_mean_field_ddt,'Ax',256,20);
% 4
% pdfA_ddt(2) = mkpdf5(trajs_conc_without_mean_field_ddt,'Ay',256,20);
% 5
% pdfA_ddt(3) = mkpdf5(trajs_conc_without_mean_field_ddt,'Az',256,20);
% 6

pdfVabs_fullg = mkpdf5(trajs_conc_with_mean_field_fullg,'Vabs',256,10);
pdfVabs_ddt = mkpdf5(trajs_conc_with_mean_field_ddt,'Vabs',256,10);

pdfAabs_fullg = mkpdf5(trajs_conc_with_mean_field_fullg,'Aabs',256,20);
pdfAabs_ddt = mkpdf5(trajs_conc_with_mean_field_ddt,'Aabs',256,20);

pdfV_fullg(1) = mkpdf5(trajs_conc_with_mean_field_fullg,'Vx',256,10);
1
pdfV_fullg(2) = mkpdf5(trajs_conc_with_mean_field_fullg,'Vy',256,10);
2
pdfV_fullg(3) = mkpdf5(trajs_conc_with_mean_field_fullg,'Vz',256,10);
3

pdfV_ddt(1) = mkpdf5(trajs_conc_with_mean_field_ddt,'Vx',256,10);
1
pdfV_ddt(2) = mkpdf5(trajs_conc_with_mean_field_ddt,'Vy',256,10);
2
pdfV_ddt(3) = mkpdf5(trajs_conc_with_mean_field_ddt,'Vz',256,10);
3

pdfA_fullg(1) = mkpdf5(trajs_conc_with_mean_field_fullg,'Ax',256,20);
4
pdfA_fullg(2) = mkpdf5(trajs_conc_with_mean_field_fullg,'Ay',256,20);
5
pdfA_fullg(3) = mkpdf5(trajs_conc_with_mean_field_fullg,'Az',256,20);
6

pdfA_ddt(1) = mkpdf5(trajs_conc_with_mean_field_ddt,'Ax',256,20);
4
pdfA_ddt(2) = mkpdf5(trajs_conc_with_mean_field_ddt,'Ay',256,20);
5
pdfA_ddt(3) = mkpdf5(trajs_conc_with_mean_field_ddt,'Az',256,20);
6

%% Get std(u)

vx_ddt = vertcat(trajs_conc_with_mean_field_ddt.Vx);
vy_ddt = vertcat(trajs_conc_with_mean_field_ddt.Vy);
vz_ddt = vertcat(trajs_conc_with_mean_field_ddt.Vz);

vx_fullg = vertcat(trajs_conc_with_mean_field_fullg.Vx);
vy_fullg = vertcat(trajs_conc_with_mean_field_fullg.Vy);
vz_fullg = vertcat(trajs_conc_with_mean_field_fullg.Vz);

std_u_fullg = std(sqrt(vx_fullg.^2 +vy_fullg.^2 +vz_fullg.^2))
std_u_ddt = std(sqrt(vx_ddt.^2 +vy_ddt.^2 +vz_ddt.^2))

ratio_std_u_fullg = std(vz_fullg)/std(vx_fullg)
ratio_std_u_ddt = std(vz_ddt)/std(vx_ddt)

%% Get moments
clc
k=3;

%disp([num2str(pdfV_fullg(k).mean) '--' num2str(pdfV_ddt(k).mean)])

%disp([num2str(pdfV_fullg(k).std) '--' num2str(pdfV_ddt(k).std)])

disp([num2str(pdfV_fullg(k).skewness) '--' num2str(pdfV_ddt(k).skewness)])

disp([num2str(pdfV_fullg(k).flatness) '--' num2str(pdfV_ddt(k).flatness)])


%%
clc
disp([num2str(pdfVabs_fullg.mean) '--' num2str(pdfVabs_ddt.mean)])
disp([num2str(pdfVabs_fullg.std) '--' num2str(pdfVabs_ddt.std)])
disp([num2str(pdfVabs_fullg.skewness) '--' num2str(pdfVabs_ddt.skewness)])
disp([num2str(pdfVabs_fullg.flatness) '--' num2str(pdfVabs_ddt.flatness)])


%% Plot Normalized PDFs
figure;

%load('trajs_with_mean_fullFOV.mat') % DONE WITH THIS DATA FOR PAPER #1

% semilogy(pdfV_fullg(1).xpdfn,pdfV_fullg(1).pdfn,'pentagram',MarkerSize=6,Color=color3(1,:), MarkerFaceAlpha = 0.5, MarkerEdgeAlpha = 0.5);hold on;
% semilogy(pdfV_fullg(2).xpdfn,pdfV_fullg(2).pdfn,'pentagram',MarkerSize=6,Color=color3(2,:),LineWidth=2);
% semilogy(pdfV_fullg(3).xpdfn,pdfV_fullg(3).pdfn,'pentagram',MarkerSize=6,Color=color3(3,:),LineWidth=2);
%semilogy(pdfV_ddt(1).xpdfn,pdfV_ddt(1).pdfn,'^',MarkerSize=6,Color=color3(1,:),LineWidth=2);hold on;
%semilogy(pdfV_ddt(2).xpdfn,pdfV_ddt(2).pdfn,'^',MarkerSize=6,Color=color3(2,:),LineWidth=2);
%semilogy(pdfV_ddt(3).xpdfn,pdfV_ddt(3).pdfn,'^',MarkerSize=6,Color=color3(3,:),LineWidth=2);


scatter(pdfV_fullg(1).xpdfn,pdfV_fullg(1).pdfn, 70, 'pentagram', 'filled', 'MarkerFaceColor', color3(1,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); hold on
scatter(pdfV_fullg(2).xpdfn,pdfV_fullg(2).pdfn, 70, 'pentagram', 'filled', 'MarkerFaceColor', color3(2,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(pdfV_fullg(3).xpdfn,pdfV_fullg(3).pdfn, 70, 'pentagram', 'filled', 'MarkerFaceColor', color3(3,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

scatter(pdfV_ddt(1).xpdfn,pdfV_ddt(1).pdfn,70, '^', 'filled', 'MarkerFaceColor', color3(1,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); hold on
scatter(pdfV_ddt(2).xpdfn,pdfV_ddt(2).pdfn,70, '^', 'filled', 'MarkerFaceColor', color3(2,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(pdfV_ddt(3).xpdfn,pdfV_ddt(3).pdfn,70, '^', 'filled', 'MarkerFaceColor', color3(3,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=3);

set(gca, 'YScale', 'log');
%set(gca,FontSize=15)
lgd = legend('$V_x$ (Earth)','$V_y$ (Earth)','$V_z$ (Earth)','$V_x$ (Microg.)','$V_y$ (Microg.)','$V_z$ (Microg.)','Gaussian','interpreter','latex',Location='south');
legendMarkers = findobj(lgd, 'type', 'patch');
set(legendMarkers, 'MarkerSize', 2500); % Set the desired marker size in the legend
%title('$PDF$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
xlabel('$\frac{V-\langle V \rangle}{std(V)}$','interpreter','latex',FontWeight='bold')
grid off
box on
axis padded
ylim([5e-7 1])
xlim([-10 10])
stop

folderout = 'pdfs';
mkdir(folderout)
savefig_FC([folderout filesep 'PDF_v'],8,6,'pdf')
savefig_FC([folderout filesep 'PDF_v'],8,6,'fig')
%%%%%%%%%

figure;

% semilogy(pdfA_fullg(1).xpdfn,pdfA_fullg(1).pdfn,'pentagram',MarkerSize=6,Color=color3(1,:),LineWidth=2); hold on
% semilogy(pdfA_fullg(2).xpdfn,pdfA_fullg(2).pdfn,'pentagram',MarkerSize=6,Color=color3(2,:),LineWidth=2);
% semilogy(pdfA_fullg(3).xpdfn,pdfA_fullg(3).pdfn,'pentagram',MarkerSize=6,Color=color3(3,:),LineWidth=2);
% 
% semilogy(pdfA_ddt(1).xpdfn,pdfA_ddt(1).pdfn,'^',MarkerSize=6,Color=color3(1,:),LineWidth=2); hold on
% semilogy(pdfA_ddt(2).xpdfn,pdfA_ddt(2).pdfn,'^',MarkerSize=6,Color=color3(2,:),LineWidth=2);
% semilogy(pdfA_ddt(3).xpdfn,pdfA_ddt(3).pdfn,'^',MarkerSize=6,Color=color3(3,:),LineWidth=2);

scatter(pdfA_fullg(1).xpdfn,pdfA_fullg(1).pdfn, 70, 'pentagram', 'filled', 'MarkerFaceColor', color3(1,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); hold on
scatter(pdfA_fullg(2).xpdfn,pdfA_fullg(2).pdfn, 70, 'pentagram', 'filled', 'MarkerFaceColor', color3(2,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(pdfA_fullg(3).xpdfn,pdfA_fullg(3).pdfn, 70, 'pentagram', 'filled', 'MarkerFaceColor', color3(3,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

scatter(pdfA_ddt(1).xpdfn,pdfA_ddt(1).pdfn,70, '^', 'filled', 'MarkerFaceColor', color3(1,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); hold on
scatter(pdfA_ddt(2).xpdfn,pdfA_ddt(2).pdfn,70, '^', 'filled', 'MarkerFaceColor', color3(2,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(pdfA_ddt(3).xpdfn,pdfA_ddt(3).pdfn,70, '^', 'filled', 'MarkerFaceColor', color3(3,:), 'MarkerEdgeColor', color3(1,:), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=3);

ylabel('$PDF(\frac{A-\langle A \rangle}{std(A)})$','interpreter','latex',FontWeight='bold')
legend('$A_x$ (Earth)','$A_y$ (Earth)','$A_z$ (Earth)','$A_x$ (Microg.)','$A_y$ (Microg.)','$A_z$ (Microg.)','Gaussian','interpreter','latex',Location='south');
xlabel('$\frac{A-\langle A \rangle}{std(A)}$','interpreter','latex',FontWeight='bold')
grid off
set(gca, 'YScale', 'log');
box on
%xlim([-5 5])
ylim([5e-7 1])
xlim([-10 10])


folderout = 'pdfs';
mkdir(folderout)
savefig_FC([folderout filesep 'PDF_a'],8,6,'pdf')
savefig_FC([folderout filesep 'PDF_a'],8,6,'fig')

%% Table with moments of distribution
maketable(pdfA_fullg,pdfV_fullg,[],[],folderout)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S2L

n=2;
start_idx = n * 5e4 + 1;
end_idx = min((n + 1) * 5e4, length(trajs_conc_without_mean_field_fullg));

S2L_fullg(1)= structFunc_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Vx',2);
S2L_fullg(2)= structFunc_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Vy',2);
S2L_fullg(3)= structFunc_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Vz',2);

j=4;
start_idx = j * 5e4 + 1;
end_idx = min((j + 1) * 5e4, length(trajs_conc_without_mean_field_ddt));

S2L_ddt(1)= structFunc_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Vx',2);
S2L_ddt(2)= structFunc_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Vy',2);
S2L_ddt(3)= structFunc_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Vz',2);

%% Plot
figure;
loglog(S2L_fullg(1).tau/Fs,S2L_fullg(1).mean,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on
loglog(S2L_fullg(2).tau/Fs,S2L_fullg(2).mean,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);
loglog(S2L_fullg(3).tau/Fs,S2L_fullg(3).mean,'d',MarkerSize=3,Color=color3(3,:),LineWidth=2);

loglog(S2L_ddt(1).tau/Fs,S2L_ddt(1).mean,'ro',MarkerSize=3,LineWidth=2);hold on
loglog(S2L_ddt(2).tau/Fs,S2L_ddt(2).mean,'go',MarkerSize=3,LineWidth=2);
loglog(S2L_ddt(3).tau/Fs,S2L_ddt(3).mean,'bo',MarkerSize=3,LineWidth=2);


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

stop
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













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eulerian 2-point statistics -- WITH MEAN FIELD

%for j=1:numel(trajs_conc); trajs_conc(j).Tf = trajs_conc(j).t_sec_abs; end % rename Tf field
%for n=1:100
n=5;
start_idx = n * 5e4 + 1;
end_idx = min((n + 1) * 5e4, length(trajs_conc_with_mean_field_fullg));
[eulerStats_fullg, pair_fullg] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_with_mean_field_fullg(start_idx:end_idx),[0.5 100],30,'off');

j=2;
start_idx = j * 5e4 + 1;
end_idx = min((j + 1) * 5e4, length(trajs_conc_with_mean_field_ddt));
[eulerStats_ddt, pair_ddt] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_with_mean_field_ddt(start_idx:end_idx),[0.5 100],15,'off');

%stop

%% Plot integral of R = 1-S2/sigma = Integral Scale
figure(999);clf

Vx = vertcat(trajs_conc_with_mean_field_fullg.Vx);
Vy = vertcat(trajs_conc_with_mean_field_fullg.Vy);
Vz = vertcat(trajs_conc_with_mean_field_fullg.Vz);

sigma2 = (std(Vx).^2+std(Vy).^2+std(Vz).^2)/3;

 % r = eulerStats_fullg.Ruur;
 % R = eulerStats_fullg.Ruu;

r = eulerStats_fullg.r;
R = 1-eulerStats_fullg.Splong{1,2}/(2*sigma2);

cumulative_integral = cumtrapz(r, R);

plot(r,cumulative_integral,'ro-')

%% Plot S2 with epsilon inset
figure(3);clf

loglog(eulerStats_fullg.r,eulerStats_fullg.Splong{1,2}./1e6,'ro-',MarkerSize=8,Color=color3(1,:),LineWidth=2);
hold on
loglog(eulerStats_ddt.r,eulerStats_ddt.Splong{1,2}./1e6,'go-',MarkerSize=8,Color=color3(2,:),LineWidth=2);

r = linspace(1,9,100);
loglog(r,7e3/1e6*r.^(2/3),'--',Color=color3(3,:),LineWidth=2)

ylabel('$S_2^{E,\parallel} (m^2 s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded

legend('Earth','Microgravity','$r^{2/3}$','interpreter','latex',Location='northwest',FontSize=20)

savefig_FC([folderout filesep 'S2_and_eps_with_meanfield'],8,6,'pdf')% just to change font

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot epsilon
insetAxes = axes('Position', [0.55, 0.2, 0.3, 0.2]); % [left, bottom, width, height]

Ckolomogrov = 2.1;
loglog(eulerStats_fullg.r,(eulerStats_fullg.Splong{1,2}./1e6./Ckolomogrov).^(1.5)./eulerStats_fullg.r','d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);
hold on
loglog(eulerStats_ddt.r,(eulerStats_ddt.Splong{1,2}./1e6./Ckolomogrov).^(1.5)./eulerStats_ddt.r','d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);

r = linspace(1,9,100);
%%% compute epsilon
epsilon_fullg = (eulerStats_fullg.Splong{1,2}./1e6./Ckolomogrov).^(1.5)./eulerStats_fullg.r';
epsilon_fullg = mean(epsilon_fullg(3:8));

epsilon_ddt = (eulerStats_ddt.Splong{1,2}/1e6./Ckolomogrov).^(1.5)./eulerStats_ddt.r';
epsilon_ddt = mean(epsilon_ddt(3:8));
%%%
loglog(r,ones(size(r)).*epsilon_fullg,'--',Color=color3(3,:),LineWidth=2);hold on
loglog(r,ones(size(r)).*epsilon_ddt,'--',Color=color3(3,:),LineWidth=2)

set(gca,FontSize=20)
%ylabel('$\epsilon = (S_2^{E,\parallel}/C_k)^{3/2}\cdot r^{-1}$','interpreter','latex',FontWeight='bold',FontSize=20)
ylabel('$\epsilon$ (m^2 s^{-3})','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r (mm)$','interpreter','latex')
%ylim([4e4 1e5])
xlim([0.6 20])
xticks([1e0 1e1 2e1])
legend('','','$\epsilon$','interpreter','latex',Location='southeast',FontSize=20)
%grid on


savefig_FC([folderout filesep 'S2_and_eps_withmeanfield'],8,6,'pdf')
savefig_FC([folderout filesep 'S2_and_eps_withmeanfield'],8,6,'fig')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tracer: find optimal ndt for Correlation function -- dt method

nmaxdt = 10;
nmaxtau = 10;
find_optimal_ndt(trajs_conc_without_mean_field_ddt(1:1e5),nmaxdt,nmaxtau,'X',1)
%find_optimal_ndt(trajs_conc_without_mean_field_ddt,nmaxdt,nmaxtau,'Y',1)
%find_optimal_ndt(trajs_conc_without_mean_field_ddt,nmaxdt,nmaxtau,'Z',1)
clear nmaxtau nmaxtau

%% Tracer: correlation function -- dt method (denosied)
ndts = [6 6 6]; % start points of correlation function ndt
ndtl = [10 10 10]; % length of correlation function ndt

disp('calculating Correlation functions -- dt method')
[tau_fullg,corrv_fullg,corra_fullg] = dtCorr(trajs_conc_without_mean_field_fullg(1:5e4),ndts,ndtl,Fs);
1
[tau_ddt,corrv_ddt,corra_ddt] = dtCorr(trajs_conc_without_mean_field_ddt(1:5e4),ndts,ndtl,Fs);
2

LagragianStats_fullg.ndts = ndts;
LagragianStats_fullg.ndtl = ndtl;
LagragianStats_fullg.tau = tau_fullg;
LagragianStats_fullg.corrv = corrv_fullg;
LagragianStats_fullg.corra = corra_fullg;

LagragianStats_ddt.ndts = ndts;
LagragianStats_ddt.ndtl = ndtl;
LagragianStats_ddt.tau = tau_ddt;
LagragianStats_ddt.corrv = corrv_ddt;
LagragianStats_ddt.corra = corra_ddt;

clear ndts ndtl

%% Tracer: plot Correlation function -- dt method (denosied)
fields = fieldnames(tau);
figure(10);clf
subplot(2,1,1)
% fullg
% for kfield = 1:numel(fields)
%     f = fields{kfield};
%     plot(tau_fullg.(f),corrv_fullg.(f)/corrv_fullg.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
% end
% ddt
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau_ddt.(f),corrv_ddt.(f)/corrv_ddt.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end

yline(0,LineWidth = 2)
%legend('$x$','$y$','$z$','Location','northeast','FontWeight','bold','FontSize',20)
%xlabel('$\tau(s)$','interpreter','latex','FontWeight','bold','FontSize',20)
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$','interpreter','latex','FontWeight','bold','FontSize',20)
grid off;

subplot(2,1,2)
% fullg
% for kfield = 1:numel(fields)
%     f = fields{kfield};
%     plot(tau_fullg.(f),corra_fullg.(f)/corra_fullg.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
% end
% ddt
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau_ddt.(f),corra_ddt.(f)/corra_ddt.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end

grid off;
yline(0,LineWidth = 2)
legend('$x$','$y$','$z$','Location','northeast','FontWeight','bold','FontSize',20,'interpreter','latex')
xlabel('$\tau(s)$','interpreter','latex','FontWeight','bold','FontSize',20);
ylabel('$\frac{\langle a(t)a(t+\tau) \rangle} {\langle a^2(t) \rangle}$','interpreter','latex','FontWeight','bold','FontSize',20);
xlim([0 0.04])

stop
folderout = 'corr';
mkdir(folderout)
savefig_FC([folderout 'corr_dt'],8,6,'pdf')
savefig_FC([folderout 'corr_dt'],8,6,'fig')

save('output_post_processing.mat','Ruu','Raa','Ruufit','Raafit','-append')





%% Other Way - Doesn't Work

%%% Velocity and Acceleration Correlations

n=2;
start_idx = n * 5e4 + 1;
end_idx = min((n + 1) * 5e4, length(trajs_conc_without_mean_field_fullg));

start_idx = 1;
end_idx = numel(trajs_conc_without_mean_field_fullg);

Ruu_fullg(1) = xcorr_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Vx',1);
Ruu_fullg(2) = xcorr_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Vy',1);
Ruu_fullg(3) = xcorr_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Vz',1);
Raa_fullg(1) = xcorr_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Ax',1);
Raa_fullg(2) = xcorr_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Ay',1);
Raa_fullg(3) = xcorr_struct(trajs_conc_without_mean_field_fullg(start_idx:end_idx),'Az',1);

j=4;
start_idx = j * 5e4 + 1;
end_idx = min((j + 1) * 5e4, length(trajs_conc_without_mean_field_ddt));

Ruu_ddt(1) = xcorr_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Vx',1);
Ruu_ddt(2) = xcorr_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Vy',1);
Ruu_ddt(3) = xcorr_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Vz',1);
Raa_ddt(1) = xcorr_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Ax',1);
Raa_ddt(2) = xcorr_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Ay',1);
Raa_ddt(3) = xcorr_struct(trajs_conc_without_mean_field_ddt(start_idx:end_idx),'Az',1);


%% Tracer: fit correlation  -- Thomas's
addpath(genpath('/Users/fcb/Documents/GitHub/Cheng_old'));

nCorrFitV = 50;
nCorrFitA = 20;

% seems 2layers works better for velocity
% while infinite layers works for accerleartion
if2layersV = 2;
if2layersA = 99;

% sometimes the fit doesn't work, set the bounded option to 0 to let it work
ifboundedV = [1 1 1];
ifboundedA = [1 1 1];

Ruufit(1) = correlationFit(Ruu_fullg(1),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(1));
Ruufit(2) = correlationFit(Ruu_fullg(2),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(2));
Ruufit(3) = correlationFit(Ruu_fullg(3),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(3));
Raafit(1) = correlationFit(Raa_fullg(1),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(1));
Raafit(2) = correlationFit(Raa_fullg(2),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(2));
Raafit(3) = correlationFit(Raa_fullg(3),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(3));
clear nCorrFitV nCorrFitA
clear if2layersV  if2layersA
clear ifboundedV ifboundedA

%% Tracer: plot Correlation function (fit)
% figure;
% main plot: zoom in
f1 = figure(11); clf
tiledlayout(2,1)
nexttile
plot(Ruu_fullg(1).tau/Fs,Ruu_fullg(1).mean/Ruu_fullg(1).mean(1),'-d',Color=color3(1,:),MarkerSize=1);hold on
plot(Ruu_fullg(2).tau/Fs,Ruu_fullg(2).mean/Ruu_fullg(2).mean(1),'-d',Color=color3(2,:),MarkerSize=1);
plot(Ruu_fullg(3).tau/Fs,Ruu_fullg(3).mean/Ruu_fullg(3).mean(1),'-d',Color=color3(3,:),MarkerSize=1);
% plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:));hold on
% plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:))
% plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:))

legend('$x$','$y$','$z$',Location='eastoutside')
title('$R_{uu}$')
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$');
xlabel('$\tau$/s')
grid on
axis tight

nexttile
plot(Raa_fullg(1).tau/Fs,Raa_fullg(1).mean/Raa_fullg(1).mean(1),'-o',Color=color3(1,:),MarkerSize=1);hold on
plot(Raa_fullg(2).tau/Fs,Raa_fullg(2).mean/Raa_fullg(2).mean(1),'-o',Color=color3(2,:),MarkerSize=1);
plot(Raa_fullg(3).tau/Fs,Raa_fullg(3).mean/Raa_fullg(3).mean(1),'-o',Color=color3(3,:),MarkerSize=1);
% plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:))
% plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:))
% plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:))

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








































%% plot S2
figure(1);clf

loglog(eulerStats_fullg.r,eulerStats_fullg.Splong{1,2},'ro-',MarkerSize=8,Color=color3(1,:),LineWidth=2);
hold on
loglog(eulerStats_ddt.r,eulerStats_ddt.Splong{1,2},'go-',MarkerSize=8,Color=color3(2,:),LineWidth=2);

r = linspace(0.4,100,100);
loglog(r,5e3*r.^(2/3),'--',Color=color1,LineWidth=2)

%legend('$(S_2^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
%title('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^{E,\parallel} (m^2 s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded


%folderout = 'S2_euler';
%mkdir(folderout)
savefig_custom([folderout filesep 'S2'],8,6,'pdf')
savefig_custom([folderout filesep 'S2'],8,6,'fig')


%% plot epsilon

Ckolomogrov = 2.1;

figure(2);clf

loglog(eulerStats_fullg.r,(eulerStats_fullg.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats_fullg.r','d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);
hold on
loglog(eulerStats_ddt.r,(eulerStats_ddt.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats_ddt.r','d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);


set(gca,FontSize=15)
legend('Earth','Microgravity','interpreter','latex',Location='best',FontSize=12)
%title('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$\epsilon = (S_2^{E,\parallel}/C_k)^{3/2}\cdot r^{-1}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded


%folderout = 'epsilon';
%mkdir(folderout)
savefig_custom([folderout filesep 'epsilon'],8,6,'pdf')
savefig_custom([folderout filesep 'epsilon'],8,6,'fig')

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
