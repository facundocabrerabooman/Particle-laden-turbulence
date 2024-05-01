%% Average and filter calibration images 
close all, clear, clc 

cam = 3;
for frame=1:6

%%%%%%%%%
pathcalib = ['I:\Calib\calib_4-30-24\raw\Camera' num2str(cam) filesep num2str(frame)];

pathout = 'I:\Calib\calib_4-30-24\Calib_Ims_4-30-24\';

fname = ['cam' num2str(cam) '_frame_preproc_00000' num2str(frame)];

mkdir(pathout)
cd(pathout) 
image_list = dir([pathcalib filesep '*.tiff']);

imsum = zeros(512,1280);

for i=5:numel(image_list)
    im = imread([pathcalib filesep image_list(i).name]);
    imsum = imadd(double(im),imsum);  
end

imaver = uint16(imsum./numel(image_list));
%imshow(imaver)

% roi_mask = roipoly();
% imaver(~roi_mask) = 0;

imaver = imaver.*(5.3e4/65500);

st=strel('disk',1);
imaver=imopen(imaver,st);
%imaver=imadjust(imaver);

imshow(imaver)

imwrite(uint16(imaver),[pathout filesep fname '.tiff'])
end



%% Average and filter calibration images _ Front Light images old
close all, clear, clc

cam = 3;
for frame=1:6

%%%%%%%%%
pathcalib = ['I:\Calib\TrCer_1000\Calib_4-18-24\raw\Camera' num2str(cam) filesep num2str(frame)];

pathout = 'I:\Calib\TrCer_1000\Calib_4-18-24\';

fname = ['cam' num2str(cam) '_frame_preproc_00000' num2str(frame)];

mkdir(pathout)
cd(pathout)
image_list = dir([pathcalib filesep '*.tiff']);


imsum = zeros(512,1280);

 %im = imread([pathcalib filesep image_list(5).name]);
 %imshow(im)
% h = imfreehand;

for i=5:numel(image_list)
    
    im = imread([pathcalib filesep image_list(i).name]);
    
    %thr = 1.5e4; % cam1
    %thr = 4e4; % cam2
    %thr = 2.8e4; % cam3

    %im(im<thr) = 0;
    %im(im>thr) = 65535;

    imsum = imadd(double(im),imsum);
    
end

imaver = uint16(imsum./numel(image_list));

st=strel('disk',1);
imaver=imopen(imaver,st);
imaver=imadjust(imaver);

imaver = imaver.*(5.3e4/65500);

imshow(imaver)

%cam1 = 65526
%cam2 = 53000
%cam3 = 53000
imwrite(uint16(imaver),[pathout filesep fname '.tiff'])
end
stop


%% Filter averaged frames

pathcalib = ['/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop ' ...
    'Tower Multiphase Flow Project/Data/bigbiginertial/calibs/calib_27nov_b/Camera1/'];

pathout = ['/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop ' ...
    'Tower Multiphase Flow Project/Data/bigbiginertial/calibs/calib_27nov_b/mean_filtered/'];
mkdir(pathout)

fname = 'cam1_frame_preproc_000001';


im = imread([pathcalib filesep fname '.tiff']);
im(im<1500)=0;

st1 = strel('disk', 10);
st2 = strel('disk', 5);

imo=imopen(im,st2);
imo=imclose(imo,st1);



%% Cut calibration frames

pathcalib = '/Users/FC/Aux_files/Calib_26oct23/';
pathout = '/Users/FC/Aux_files/Calib_26oct23/cut/';
mkdir(pathout)
image_list = dir([pathcalib filesep '*.tiff']);

for i=1:numel(image_list)
    
    im = imread([pathcalib filesep image_list(i).name]);
    im = im(285:796,321:1600); 
    imwrite(uint16(im),[pathout filesep image_list(i).name])


end

