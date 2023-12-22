tic
close all
clear 
clc

% input
fpath0 = 'C:\Chronos_Footages\chain_1107_3kHz_anagain6dB';
Istart = 2; % better to ignore the 1st image


fpath = [fpath0 '\raw'];
% get the folder contents
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));

%%
ifsave = 1;

stdil = strel('disk',7);
%%
if ifsave == 1
    preproc_dirp = [fpath0 '\preproc\'];
    mkdir(preproc_dirp)
end

image_list = dir([fpath '\*.tiff']);
img_num = size(image_list,1)/3;

for k=Istart:img_num
% for k=100:100
        
    %% cam1
    fname=[fpath '\' 'cam1_frame_' num2str(k,'%06d') '.tiff'];
    Im1 = imread(fname); 
    Im2  = imopen(Im1,stdil);
    [C, R] = imfindcircles(Im2, [15 20],'ObjectPolarity','dark','Sensitivity',0.95);
%     figure;imagesc(Im1);hold on;viscircles(C,R)
    mask = zeros(size(Im1,1),size(Im1,2));
    for i = 1:size(C,1)
        roi=images.roi.Circle('Center',C(i,:),'Radius',R(i));
        mask = mask + createMask(roi,size(Im1,1),size(Im1,2));
    end
    Imp1 = immultiply(Im1,cast(mask,class(Im1)));
    Imp1 = imgaussfilt(Imp1,1);

    %% cam2
    fname=[fpath '\' 'cam2_frame_' num2str(k,'%06d') '.tiff'];
    Im1 = imread(fname); 
    Im2  = imopen(Im1,stdil);
    [C, R] = imfindcircles(Im2, [15 20],'ObjectPolarity','dark','Sensitivity',0.95);
%     figure;imagesc(Im1);hold on;viscircles(C,R)
    mask = zeros(size(Im1,1),size(Im1,2));
    for i = 1:size(C,1)
        roi=images.roi.Circle('Center',C(i,:),'Radius',R(i));
        mask = mask + createMask(roi,size(Im1,1),size(Im1,2));
    end
    Imp2 = immultiply(Im1,cast(mask,class(Im1)));
    Imp2 = imgaussfilt(Imp2,1);

    %% cam3
    fname=[fpath '\' 'cam3_frame_' num2str(k,'%06d') '.tiff'];
    Im1 = imread(fname); 
    Im2  = imopen(Im1,stdil);
    [C, R] = imfindcircles(Im2, [15 20],'ObjectPolarity','dark','Sensitivity',0.95);
%     figure;imagesc(Im1);hold on;viscircles(C,R)
    mask = zeros(size(Im1,1),size(Im1,2));
    for i = 1:size(C,1)
        roi=images.roi.Circle('Center',C(i,:),'Radius',R(i));
        mask = mask + createMask(roi,size(Im1,1),size(Im1,2));
    end
    Imp3 = immultiply(Im1,cast(mask,class(Im1)));
    Imp3 = imgaussfilt(Imp3,1);

    fnameo1 = ['cam1_frame_preproc_' num2str(k,'%06d') '.tiff'];
    fnameo2 = ['cam2_frame_preproc_' num2str(k,'%06d') '.tiff'];
    fnameo3 = ['cam3_frame_preproc_' num2str(k,'%06d') '.tiff'];
    
    if ifsave == 1
        imwrite(uint16(Imp1),[preproc_dirp '\' fnameo1]);
        imwrite(uint16(Imp2),[preproc_dirp '\' fnameo2]);
        imwrite(uint16(Imp3),[preproc_dirp '\' fnameo3]);
    end
end
toc