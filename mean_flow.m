%%
clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Volumes/landau1/Tracers/ddt_filtered_w10';
folderout = folderin;
cd(folderin)

Fs=2990; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

%%
load('traj_conc_tracers_ddt.mat')
tracklong_tracer = trajs_conc; clear trajs_conc
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
