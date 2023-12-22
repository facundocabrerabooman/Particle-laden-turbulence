% check if the tracers are sufficient
nexp = '01';

ind = 300;

fin = 'D:\Chronos_Footage\1.1g';
im1 = imread([fin '\' nexp '\cam1_frame_' numstr(ind) '.tiff']);
im2 = imread([fin '\' nexp '\cam2_frame_' numstr(ind) '.tiff']);
im3 = imread([fin '\' nexp '\cam3_frame_' numstr(ind) '.tiff']);

figure;imagesc(im1);
figure;imagesc(im2);
figure;imagesc(im3);