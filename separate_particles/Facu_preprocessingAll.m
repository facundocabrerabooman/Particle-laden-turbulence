close all
clear
clc
%pause(60*60)
% input
fpath0 = 'I:\Calib\calib_4-30-24\bench_test\';
cd(fpath0)

part_radius = 6;
 
fpath = [fpath0 filesep 'raw'];
% get the folder contents
addpath(genpath('C:\Users\meatlab1-admin\Documents\facundo\Particle-laden-turbulence'));
%dfolders = FunSubfolder(fpath);
dfolders = [];

%%% preproc & split
datasets = {'fullg','ddt','dec'};

for j = 1%:3
dataset = datasets{j};

    %image_list = dir([subfolders filesep '*.tiff']);
    image_list = dir([fpath filesep 'Camera1' filesep '*.tiff']);
    img_num = size(image_list,1);
 
    %Istart = 2;
    %Iend = img_num-10;
if pi==pi
    if dataset(1) == 'f'
        Istart = 2;
        Iend = img_num - round(2990*4.1);
    elseif dataset == 'ddt'
        Istart = img_num - round(2990*4.1);
        Iend = img_num - round(2990*2);
    else
        Istart = img_num - round(2990*2);
        Iend = img_num;
    end
end
      Istart = 2;
      Iend=1000;
      disp('careful here')
     % pause


    if numel(dfolders)~=0
        subfolders = [fpath filesep dfolders(i).name];
        preproc_dirt = [fpath0 filesep 'pre_proc' filesep 'preproc_tracer_' dataset filesep dfolders(i).name ];
        preproc_dirp = [fpath0 filesep 'pre_proc' filesep 'preproc_particle_' dataset filesep dfolders(i).name ];
    else
        subfolders = fpath;
        preproc_dirt = [fpath0 filesep 'pre_proc' filesep 'preproc_tracer_' dataset];
        preproc_dirp = [fpath0 filesep 'pre_proc' filesep 'preproc_particle_' dataset];
    end
    mkdir(preproc_dirt)
    mkdir(preproc_dirp)


    %%% get Backgrounds
    %cam 1
    bkg1_originalSize = getBkg(subfolders,'Camera1',Istart,Iend,1,[]);
    
    %bkg1_originalSize = imread('/Users/FC/Downloads/im0_c1.tiff');
    %     bkg1 = cast(zeros(1080,1920),class(bkg1_originalSize));
    %     bkg1(285:796,321:1600) = bkg1_originalSize;
    %cam2
    bkg2_originalSize = getBkg(subfolders,'Camera2',Istart,Iend,1,[]);
    
    %bkg2_originalSize = imread('/Users/FC/Downloads/im0_c2.tiff');
    %     bkg2 = cast(zeros(1080,1920),class(bkg2_originalSize));
    %     bkg2(285:796,321:1600) = bkg2_originalSize;

    %cam3
    bkg3_originalSize = getBkg(subfolders,'Camera3',Istart,Iend,1,[]);
    %disp('no cam3')
    %bkg3_originalSize = imread('/Users/FC/Downloads/im0_c3.tiff');
    %     bkg3 = cast(zeros(1080,1920),class(bkg3_originalSize));
    %     bkg3(285:796,321:1600) = bkg3_originalSize;
    counter = 0;

    for k=Istart:Iend
        %disp('CAREFUL')
        counter = counter+1;
        k/Iend

        %%% cam1
        fname=[subfolders filesep 'Camera1' filesep 'frame_' num2str(k,'%06d') '.tiff'];
        Im1_originalSize = imread(fname);

        intensity_thr = 1e3;
        [~,Im1t,Im1p]=Facu_preprocessing(Im1_originalSize,intensity_thr,part_radius,1,bkg1_originalSize);
        fnameo = ['cam1_frame_preproc_' num2str(counter,'%06d') '.tiff'];
        imwrite(uint16(Im1t),[preproc_dirt filesep fnameo])
        imwrite(uint16(Im1p),[preproc_dirp filesep fnameo]);
%                 figure(10);
%                 subplot(2,1,1);imagesc(Im1_originalSize);axis equal
%                 subplot(2,1,2);imagesc(Im1p);axis equal
%                 pause(0.05)


        %%% cam2
        fname=[subfolders filesep 'Camera2' filesep 'frame_'  num2str(k,'%06d') '.tiff'];
        Im2_originalSize = imread(fname);
        %         Im2 = cast(zeros(1080,1920),class(Im2_originalSize));
        %         Im2(285:796,321:1600) = Im2_originalSize;
        intensity_thr = 1e3;
        [~,Im2t,Im2p]=Facu_preprocessing(Im2_originalSize,intensity_thr,part_radius,2,bkg2_originalSize);
        fnameo = ['cam2_frame_preproc_' num2str(counter,'%06d') '.tiff'];
       imwrite(uint16(Im2t),[preproc_dirt filesep fnameo]);
        imwrite(uint16(Im2p),[preproc_dirp filesep fnameo]);
%                figure(10);
%                  subplot(2,1,1);imagesc(Im2_originalSize);axis equal
%                  subplot(2,1,2);imagesc(Im2p);axis equal
%                  pause(0.05)

if pi==pi
        %%% cam3
        fname=[subfolders filesep 'Camera3' filesep 'frame_' num2str(k,'%06d') '.tiff'];
        Im3_originalSize = imread(fname);
        %         Im3 = cast(zeros(1080,1920),class(Im3_originalSize));
        %         Im3(285:796,321:1600) = Im3_originalSize;
        intensity_thr = 1e3;
        [~,Im3t,Im3p] =Facu_preprocessing(Im3_originalSize,intensity_thr,part_radius,3,bkg3_originalSize);
        fnameo = ['cam3_frame_preproc_' num2str(counter,'%06d') '.tiff'];
        imwrite(uint16(Im3t),[preproc_dirt filesep fnameo]);
        imwrite(uint16(Im3p),[preproc_dirp filesep fnameo]);
%                          figure(10);
%                  subplot(2,1,1);imagesc(Im3_originalSize);axis equal
%                  subplot(2,1,2);imagesc(Im3p);axis equal
%                  pause(0.05)
end
    end
end
stop
%% Test camera 1

if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im1_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im1p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r')


    subplot(2,1,2)
    imagesc(Im1p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r')

    figure;imagesc(Im1t);axis equal;hold on
    %viscircles(C,ones(size(C,1),1).*5,'Color','r')

end
%% Test camera 2
if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im2_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im2p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r');


    subplot(2,1,2)
    imagesc(Im2p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r');

    figure;imagesc(Im2t);axis equal;hold on

end
%% Test camera 3
if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im3_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im3p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r');


    subplot(2,1,2)
    imagesc(Im3p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r');

    figure;imagesc(Im3t);axis equal;hold on
end