clc;clear;close all;
fin = 'F:\Chronos_Footage\chain_0927_3kHz_anagain12dB\raw';
fout = 'F:\Chronos_Footage\chain_0927_3kHz_anagain12dB';


im1 = imread([fin filesep 'cam1_frame_000001.tiff']);
[C1,R1] = imfindcircles(im1,[15 25],"ObjectPolarity","dark");
figure;imagesc(im1);hold on;viscircles(C1,R1);
C1 = sortrows(C1);
plot(C1(:,1),C1(:,2),'g+-')
p1 = polyfit(C1(:,1),C1(:,2),1);
plot([min(xlim):1:max(xlim)],p1(1)*[min(xlim):1:max(xlim)]+p1(2),'k.')


im2 = imread([fin filesep 'cam2_frame_000001.tiff']);
[C2,R2] = imfindcircles(im2,[15 25],"ObjectPolarity","dark");
figure;imagesc(im2);hold on;viscircles(C2,R2);
C2 = sortrows(C2);
plot(C2(:,1),C2(:,2),'g+-')
p2 = polyfit(C2(:,1),C2(:,2),1);
plot([min(xlim):1:max(xlim)],p2(1)*[min(xlim):1:max(xlim)]+p2(2),'k.')

im3 = imread([fin filesep 'cam3_frame_000001.tiff']);
[C3,R3] = imfindcircles(im3,[15 25],"ObjectPolarity","dark");
figure;imagesc(im3);hold on;viscircles(C3,R3);
C3 = sortrows(C3);
plot(C3(:,1),C3(:,2),'g+-')
p3 = polyfit(C3(:,1),C3(:,2),1);
plot([min(xlim):1:max(xlim)],p3(1)*[min(xlim):1:max(xlim)]+p3(2),'k.')

save([fout filesep 'vertical.mat'],'p1','p2','p3')