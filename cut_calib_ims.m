%% Cut calibration frames

pathcalib = 'C:\Users\meatlab1-admin\Downloads\New folder\4\raw\';
pathout = 'C:\Users\meatlab1-admin\Downloads\New folder\4\cut\';
mkdir(pathout)
image_list = dir([pathcalib filesep '*.tiff']);

for i=1:numel(image_list)
    
    im = imread([pathcalib filesep image_list(i).name]);
    im = im(285:796,321:1600); 
    imwrite(uint16(im),[pathout filesep image_list(i).name])


end

