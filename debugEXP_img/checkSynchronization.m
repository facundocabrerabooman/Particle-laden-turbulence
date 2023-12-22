clc;clear

fin = 'F:\Chronos_Footage\0.7g_1003_00V_3kHz_12dB\raw';

string1 = '001';

fin = [fin filesep string1];


image_list = dir([fin '\*.tiff']);
img_num = size(image_list,1)/3;

for k = 1:img_num
    Iinfo1 = imfinfo([fin '\cam1_frame_' num2str(k,'%06d') '.tiff']);
    Iinfo2 = imfinfo([fin '\cam2_frame_' num2str(k,'%06d') '.tiff']);
    Iinfo3 = imfinfo([fin '\cam3_frame_' num2str(k,'%06d') '.tiff']);

    Its.min1 = Iinfo1.DigitalCamera.DateTimeDigitized([end-4:end-3]);
    Its.sec1 = Iinfo1.DigitalCamera.DateTimeDigitized([end-1:end]);
    Its.min2 = Iinfo2.DigitalCamera.DateTimeDigitized([end-4:end-3]);
    Its.sec2 = Iinfo2.DigitalCamera.DateTimeDigitized([end-1:end]);
    Its.min3 = Iinfo3.DigitalCamera.DateTimeDigitized([end-4:end-3]);
    Its.sec3 = Iinfo3.DigitalCamera.DateTimeDigitized([end-1:end]);

    Its.ind1(k) = str2double(Its.min1)*60+str2double(Its.sec1);
    Its.ind2(k) = str2double(Its.min2)*60+str2double(Its.sec2);
    Its.ind3(k) = str2double(Its.min3)*60+str2double(Its.sec3);
end

Its.diff1 = Its.ind1-Its.ind1;
Its.diff2 = Its.ind2-Its.ind1;
Its.diff3 = Its.ind3-Its.ind1;
%%
figure;
plot(Its.ind1,'r')
hold on
plot(Its.ind2,'g')
plot(Its.ind3,'b')

figure;
plot(Its.diff1-Its.diff1(1),'r')
hold on
plot(Its.diff2-Its.diff2(1),'g')
plot(Its.diff3-Its.diff3(1),'b')
%% max/min intensity
% parfor k = 2:img_num
%     disp(num2str(floor((k/img_num)*10)/10))
%     maxIntensity1(fin,k);
% end
% Imax = maxIntensity2(fin,img_num);
% 
% %%
% figure;
% % subplot(2,1,1)
% plot(Imax.cam1,'r');hold on;plot(Imax.cam2,'g');plot(Imax.cam3,'b')
% legend('cam1','cam2','cam3')
% 
% %%
% % subplot(2,1,2);plot(double(Imax.cam1)-double(Imax.cam2),'r');
% % hold on;plot(double(Imax.cam1)-double(Imax.cam3),'g');legend('cam1-2','cam1-3')
% 
% 
% %% function
% function maxIntensity1(fin,k)
%     im1 = imread([fin '\cam1_frame_' num2str(k,'%06d') '.tiff']);
%     im2 = rgb2gray(demosaic(imread([fin '\cam2_frame_' num2str(k,'%06d') '.tiff']),'grbg'));
%     im3 = imread([fin '\cam3_frame_' num2str(k,'%06d') '.tiff']);
%     Imax.cam1 = max(max(im1));
%     Imax.cam2 = max(max(im2));
%     Imax.cam3 = max(max(im3));
%     fname = [fin '\Imax' num2str(k) '.mat'];
%     save(fname, 'Imax')
% end
% 
% function Imax = maxIntensity2(fin,img_num)
%     Imax.cam1 = [];Imax.cam2 = [];Imax.cam3 = [];
%     for i = 2:img_num
%         a = load([fin '\Imax' num2str(i) '.mat']);
%         Imax.cam1 = [Imax.cam1;a.Imax.cam1];
%         Imax.cam2 = [Imax.cam2;a.Imax.cam2];
%         Imax.cam3 = [Imax.cam3;a.Imax.cam3];
%     end
%     delete([fin,'\*.mat'])
% end