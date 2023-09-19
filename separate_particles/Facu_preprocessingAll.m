close all
clear 
clc

% input
fpath0 = '/Volumes/landau1/Experiment_MP/7p8&tracers_sept23/6sept_a';
Istart = 2; % better to ignore the 1st image

%%
fpath = [fpath0 filesep 'raw'];
% get the folder contents
addpath(genpath('/Users/FC/Documents/GitHub/Particle-laden-turbulence'));
%dfolders = FunSubfolder(fpath);
dfolders = [];

% parfor i = 1:size(dfolders,1)
for i = 1:1
    if numel(dfolders)~=0
        subfolders = [fpath filesep dfolders(i).name];
        preproc_dirt = [fpath0 filesep 'preproc' filesep 'preproc_tracer' filesep dfolders(i).name ];
        preproc_dirp = [fpath0 filesep 'preproc' filesep 'preproc_particle' filesep dfolders(i).name ];
    else
        subfolders = fpath;
        preproc_dirt = [fpath0 filesep 'preproc' filesep 'preproc_tracer' ];
        preproc_dirp = [fpath0 filesep 'preproc' filesep 'preproc_particle' ];
    end
     mkdir(preproc_dirt)
     mkdir(preproc_dirp)

    %image_list = dir([subfolders filesep '*.tiff']);
    image_list = dir([fpath filesep 'Camera1' filesep '*.tiff']);
    img_num = size(image_list,1)/3;


    %% get Backgrounds
    %cam 1
    bkg1_originalSize = getBkg(subfolders,'cam1_frame_',Istart,img_num,100,[]);
%     bkg1 = cast(zeros(1080,1920),class(bkg1_originalSize));
%     bkg1(285:796,321:1600) = bkg1_originalSize;
    %cam2
    bkg2_originalSize = getBkg(subfolders,'cam2_frame_',Istart,img_num,100,[]);
%     bkg2 = cast(zeros(1080,1920),class(bkg2_originalSize));
%     bkg2(285:796,321:1600) = bkg2_originalSize;

    %cam3
    bkg3_originalSize = getBkg(subfolders,'cam3_frame_',Istart,img_num,100,[]);
%     bkg3 = cast(zeros(1080,1920),class(bkg3_originalSize));
%     bkg3(285:796,321:1600) = bkg3_originalSize;

%     for k=Istart:img_num
    for k=300:300

        %% cam1
        fname=[subfolders filesep 'cam1_frame_' num2str(k,'%06d') '.tiff'];
        Im1_originalSize = imread(fname); 
%         Im1 = cast(zeros(1080,1920),class(Im1_originalSize));
%         Im1(285:796,321:1600) = Im1_originalSize;
        [Im01,Im1t,Im1p]=Facu_preprocessing(Im1_originalSize,1000,5,1,bkg1_originalSize);
        fnameo = ['cam1_frame_preproc_' num2str(k,'%06d') '.tiff'];
%         imwrite(uint16(Im1t),[preproc_dirt filesep fnameo]);

    
        %% cam2
        fname=[subfolders filesep 'cam2_frame_' num2str(k,'%06d') '.tiff'];
        Im2_originalSize = imread(fname);
%         Im2 = cast(zeros(1080,1920),class(Im2_originalSize));
%         Im2(285:796,321:1600) = Im2_originalSize;
        [Im02,Im2t,Im2p]=Facu_preprocessing(Im2_originalSize,1000,5,1,bkg2_originalSize);
        fnameo = ['cam2_frame_preproc_' num2str(k,'%06d') '.tiff'];
%         imwrite(uint16(Im2t),[preproc_dirt filesep fnameo]);
    
        %% cam3
        fname=[subfolders filesep 'cam3_frame_' num2str(k,'%06d') '.tiff'];
        Im3_originalSize = imread(fname);
%         Im3 = cast(zeros(1080,1920),class(Im3_originalSize));
%         Im3(285:796,321:1600) = Im3_originalSize;
        [Im03,Im3t,Im3p] =Facu_preprocessing(Im3_originalSize,1000,5,1,bkg3_originalSize);
        fnameo = ['cam3_frame_preproc_' num2str(k,'%06d') '.tiff'];
%         imwrite(uint16(Im3t),[preproc_dirt filesep fnameo]);
    end
end


%%
figure;imagesc(Im1_originalSize);axis equal;
figure;imagesc(Im01);axis equal;
figure;imagesc(Im1p);axis equal;
figure;imagesc(Im1t);axis equal;
%%
% figure;imagesc(Im2_originalSize);axis equal;
figure;imagesc(Im02);axis equal;
figure;imagesc(Im2p);axis equal;
figure;imagesc(Im2t);axis equal;
%%
% figure;imagesc(Im3_originalSize);axis equal;
figure;imagesc(Im03);axis equal;
figure;imagesc(Im3p);axis equal;
figure;imagesc(Im3t);axis equal;