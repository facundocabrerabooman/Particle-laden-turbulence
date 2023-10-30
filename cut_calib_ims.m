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

