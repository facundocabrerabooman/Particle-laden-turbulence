clear;close all;clc
fin = 'F:\Chronos_Footage\0929_1.0g_00V_3kHz_12dB';

addpath 'C:\Users\Gsu\Desktop\SDT_EXP'
dfolders = FunSubfolder(fin);


for i = 1:1:size(dfolders,1)
    subfolders = [fin '\' dfolders(i).name];
    
    im = imread([subfolders '\frame_000300.tiff']);
    figure
    imagesc(im);
    disp(dfolders(i).name)
end