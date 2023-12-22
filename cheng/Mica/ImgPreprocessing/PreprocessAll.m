tic
close all
clear 
clc

set(groot,'DefaultFigureWindowStyle','docked')

% input
fpath0 = 'F:\Chronos_Footage\1.2g_1031_15V_3kHz_6dB';
Istart = 2; % better to ignore the 1st image


fpath = [fpath0 '\raw'];
% get the folder contents
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));
dfolders = FunSubfolder(fpath);

%%
ifparticle = 1;
iftracer  = 1;
iffilter = 0; % if not use median filter and gaussian filter, set to zero -> only susbtract bkg and separate treacer from intertial particle 
ifsave = 1;

part_radius = 7;

% active when iffilter == 1
thres_noise1 = 2000;
thres_noise2 = 4000;
thres_noise3 = 4000;

%
thres_part1 = 2000;
thres_part2 = 2000;
thres_part3 = 2000;

if ifparticle == 0
    part_radius = 10;   
end

if ifsave == 0
    Istart = 300;
end

%%
parfor i = 1:size(dfolders,1)
% for i = 16:16
    subfolders = [fpath '\' dfolders(i).name];
    if ifsave == 1
        if ifparticle == 1
            preproc_dirp = [fpath0 '\preproc\preproc_particle\' dfolders(i).name ];
            mkdir(preproc_dirp)
        end
        if iftracer == 1
            preproc_dirt = [fpath0 '\preproc\preproc_tracer\' dfolders(i).name ];
            mkdir(preproc_dirt)
        end
    end
    
    image_list = dir([subfolders '\*.tiff']);
    img_num = size(image_list,1)/3;

    % get Backgrounds
    %cam 1
    bkg1 = getBkg(subfolders,'cam1_frame_',Istart,img_num,5,[]);
    meanI1 = double(mean(mean(bkg1(:,1:1200))));
    %cam2
    bkg2 = getBkg(subfolders,'cam2_frame_',Istart,img_num,5,[]);
    meanI2 = double(mean(mean(bkg2(:,1:1200))));
    %cam3
    bkg3 = getBkg(subfolders,'cam3_frame_',Istart,img_num,5,[]);
    meanI3 = double(mean(mean(bkg3(:,1:1200))));

    meanI = meanI3;

    if ifsave == 0
        img_num = Istart;
    end
    for k=Istart:img_num
        
        %% cam1
        fname=[subfolders '\' 'cam1_frame_' num2str(k,'%06d') '.tiff'];
        Im1 = imread(fname); 
        [Im1_sub,Im1t,Im1p]=preprocessing(Im1,thres_noise1,part_radius,thres_part1,1,bkg1,iffilter);


        %% cam2
        fname=[subfolders '\' 'cam2_frame_' num2str(k,'%06d') '.tiff'];
        Im2 = imread(fname);
        [Im2_sub,Im2t,Im2p]=preprocessing(Im2,thres_noise2,part_radius,thres_part2,1,bkg2,iffilter);
    
        %% cam3
        fname=[subfolders '\' 'cam3_frame_' num2str(k,'%06d') '.tiff'];
        Im3 = imread(fname); 
        [Im3_sub,Im3t,Im3p]=preprocessing(Im3,thres_noise3,part_radius,thres_part3,1,bkg3,iffilter);

        %%
        
        %%
        fnameo1 = ['cam1_frame_preproc_' num2str(k,'%06d') '.tiff'];
        fnameo2 = ['cam2_frame_preproc_' num2str(k,'%06d') '.tiff'];
        fnameo3 = ['cam3_frame_preproc_' num2str(k,'%06d') '.tiff'];
        if ifsave == 1
            Im1t = uint16(double(Im1t)/meanI1*meanI);
            Im1p = uint16(double(Im1p)/meanI1*meanI);
            Im2t = uint16(double(Im2t)/meanI2*meanI);
            Im2p = uint16(double(Im2p)/meanI2*meanI);
            Im3t = uint16(double(Im3t)/meanI3*meanI);
            Im3p = uint16(double(Im3p)/meanI3*meanI);
            
            if ifparticle == 1
                imwrite(uint16(Im1p),[preproc_dirp '\' fnameo1]);
                imwrite(uint16(Im2p),[preproc_dirp '\' fnameo2]);
                imwrite(uint16(Im3p),[preproc_dirp '\' fnameo3]);
            end
            if iftracer == 1
                imwrite(uint16(Im1t),[preproc_dirt '\' fnameo1]);
                imwrite(uint16(Im2t),[preproc_dirt '\' fnameo2]);
                imwrite(uint16(Im3t),[preproc_dirt '\' fnameo3]);
            end
        end
        
    end
    disp(['Folder has been Processed! ' num2str(i)])
end
toc

%%
if ifsave ==0
    close all
    figure;imagesc(bkg1);axis equal;title('bkg1');
    figure;imagesc(Im1);axis equal;title('Raw1');
    figure;imagesc(Im1_sub );axis equal;title('Cam1-Substract');
    figure;imagesc(Im1t);axis equal;title('Cam1-Tracers');
    figure;imagesc(Im1p);axis equal;title('Cam1-Particle');
    %%
    close all
    figure;imagesc(bkg2);axis equal;title('bkg2');
    figure;imagesc(Im2);axis equal;title('Raw2');
    figure;imagesc(Im2_sub);axis equal;title('Cam2-Substract');
    figure;imagesc(Im2t);axis equal;title('Cam2-Tracers');
    figure;imagesc(Im2p);axis equal;title('Cam2-Particle');
    %%
    close all
    figure;imagesc(bkg3);axis equal;title('bkg3');
    figure;imagesc(Im3);axis equal;title('Raw3');
    figure;imagesc(Im3_sub);axis equal;title('Cam3-Substract');
    figure;imagesc(Im3t);axis equal;title('Cam3-Tracers');
    figure;imagesc(Im3p);axis equal;title('Cam3-Particle');
end