clear;
clc;
close all

fin = ['F:\Chronos_Footage\1.2g_1009_00V_3kHz_12dB'];

addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'))
dfolders = FunSubfolder(fin);


for i = 1:1:size(dfolders,1)
    subfolders = [fin '\' dfolders(i).name];

    image_list1 = dir([subfolders '\*.tiff']);
    img_num1(i) = size(image_list1,1);
    disp([i img_num1(i)])
end
plot(img_num1)
mean(img_num1)/3