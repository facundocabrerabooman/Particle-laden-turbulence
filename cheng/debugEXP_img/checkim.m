clear;
clc;
close all
set(groot,'DefaultFigureWindowStyle','docked')
fin = 'C:\Chronos_Footages\constB_1.0g_1108_12V_3kHz_12dB\raw\';

for i = 1:5:100
    
    string0 = num2str(i,'%03d');
    fin0 = [fin string0];
    
    image_list1 = dir([fin0 '_cam1\*.tiff']);
    img_num1 = size(image_list1,1);
    image_list2 = dir([fin0 '_cam2\*.tiff']);
    img_num2 = size(image_list2,1);
    image_list3 = dir([fin0 '_cam3\*.tiff']);
    img_num3 = size(image_list3,1);
    
    if img_num1==img_num2 && img_num1 == img_num3
        disp(img_num1)
    else
        disp([img_num1 img_num2 img_num3])
        warning(['NOT equal'])
    end
    % fname=[fin '_cam1\frame_000300.tiff'];
    fname=[fin0 '_cam1\frame_000300.tiff'];
    im1 = imread(fname);
    subplot(3,1,1);imagesc(im1);colorbar
    title(string0)
    
    fname=[fin0 '_cam2\frame_000300.tiff'];
    im2 = imread(fname);
    subplot(3,1,2);imagesc(im2);colorbar
    title(string0)
    
    fname=[fin0 '_cam3\frame_000300.tiff'];
    im3 = imread(fname);
    subplot(3,1,3);imagesc(im3);colorbar
    title(string0)
    close all;clc
end