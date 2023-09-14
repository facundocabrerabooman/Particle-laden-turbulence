clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/FC/Documents/GitHub/Particle-laden-turbulence'));

% Set as current directory the folder with the data
folderin = '/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7a';
cd(folderin)

Fs=2990; % Frame rate

%% Set cool colors for plots
mycolormap = mycolor('#063970','#eeeee4','#e28743');
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Load data
fname = 'july7a_tracers';
load('output_post_processing','tracklong');

%% Get Velocity Fields

folderout = [folderin filesep 'fields' filesep 'outputs' filesep];
mkdir(folderout)

%%%%% Mean
[mXdt, mBdt, bins] = track2meanDxDt3DProfile(tracklong,'Xf',[4 6 8 10],[10 10 10],1,1,'x','cart');
[mYdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Yf',[4 6 8 10],[10 10 10],1,1,'y','cart');
[mZdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Zf',[4 6 8 10],[10 10 10],1,1,'z','cart');

[X,Y,Z]=meshgrid(bins{1},bins{2},bins{3});

save([folderout 'output_Vel_meanfields'],'mXdt','mYdt','mZdt','mBdt','bins','X','Y','Z')
clearvars -except Fs folderout tracklong fname folderin color3 color1
%%%%% RMS
Fs = 2990;

[mXdt mBdt bins] = track2meanDxDt3DProfile(tracklong,'Xf',[4 6 8 10],[10 10 10],1,2,'x','cart');
[mYdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Yf',[4 6 8 10],[10 10 10],1,2,'y','cart');
[mZdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Zf',[4 6 8 10],[10 10 10],1,2,'z','cart');

[X,Y,Z]=meshgrid(bins{1},bins{2},bins{3});

save([folderout 'output_Vel_rmsfields'],'mXdt','mYdt','mZdt','mBdt','bins','X','Y','Z')
clearvars -except Fs folderout tracklong fname folderin color3 color1
%% Get Acceleration Fields

folderout = [folderin filesep 'fields' filesep 'outputs' filesep];


%%%%% Mean
[mXdt, mBdt, bins] = track2meanDxDt3DProfile(tracklong,'Xf',[4 6 8 10],[10 10 10],2,1,'x','cart');
[mYdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Yf',[4 6 8 10],[10 10 10],2,1,'y','cart');
[mZdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Zf',[4 6 8 10],[10 10 10],2,1,'z','cart');

[X,Y,Z]=meshgrid(bins{1},bins{2},bins{3});

save([folderout 'output_Acc_meanfields'],'mXdt','mYdt','mZdt','mBdt','bins','X','Y','Z')
clearvars -except Fs folderout tracklong fname folderin color3 color1

%%%%% RMS

[mXdt mBdt bins] = track2meanDxDt3DProfile(tracklong,'Xf',[4 6 8 10],[10 10 10],2,2,'x','cart');
[mYdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Yf',[4 6 8 10],[10 10 10],2,2,'y','cart');
[mZdt, ~, ~] = track2meanDxDt3DProfile(tracklong,'Zf',[4 6 8 10],[10 10 10],2,2,'z','cart');

[X,Y,Z]=meshgrid(bins{1},bins{2},bins{3});

save([folderout 'output_Acc_rmsfields'],'mXdt','mYdt','mZdt','mBdt','bins','X','Y','Z')
clearvars -except Fs folderout tracklong fname folderin color3 color1

%% Plot Mean Velocity
close all, clear all, clc

pathin = '/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7a/fields/';
folderin = [pathin filesep 'outputs' filesep];
folderout = [pathin filesep 'meanvel_figures' filesep];
mkdir(folderout)

outputFileName = 'output_Vel_meanfields.mat';

plot_fields(folderin, folderout, outputFileName)

%% Plot RMS Velocity
close all, clear all, clc

pathin = '/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7a/fields/';
folderin = [pathin filesep 'outputs' filesep];
folderout = [pathin filesep 'rmsvel_figures' filesep];
mkdir(folderout)

outputFileName = 'output_Vel_rmsfields.mat';

plot_fields(folderin, folderout, outputFileName)

%% Plot Mean Acceleration
close all, clear all, clc

pathin = '/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7a/fields/';
folderin = [pathin filesep 'outputs' filesep];
folderout = [pathin filesep 'meanacc_figures' filesep];
mkdir(folderout)

outputFileName = 'output_Acc_meanfields.mat';

plot_fields(folderin, folderout, outputFileName)

%% Plot RMS Acceleration
close all, clear all, clc

pathin = '/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7a/fields/';
folderin = [pathin filesep 'outputs' filesep];
folderout = [pathin filesep 'rmsacc_figures' filesep];
mkdir(folderout)

outputFileName = 'output_Acc_rmsfields.mat';

plot_fields(folderin, folderout, outputFileName)