function [Im_sub,Imt,Imp]=preprocessing(Im,th,part_radius,gain,bkg)

% 03/2019 - Thomas Basset
%Im = imcomplement(Im);


%se=strel('disk',strel_size); %opening to remove big elements
%imo=imopen(Im,se);

Im=imsubtract(bkg,Im);
Im_sub = Im;
%Im = bpass2(Im,lnoise,part_size);
%Im = imgaussfilt(Im,part_size);
%Im=imnlmfilt(Im);

% Th = th*mean(mean(Im));
Im(Im<th)=0; %thresholding to remove the background
Im=medfilt2(Im,[3,3]);
%Im = imadjust(Im,stretchlim(Im),[]);

%Im=imsharpen(Im,'Radius',2);
%Im(Im<th)=0; %thresholding to remove the background

%se=strel('disk',1); %opening to remove little elements
%Im=imerode(Im,se);

%Imbw=imregionalmax(Im,6);
%Im = Im/mean(Im(Imbw))*gain;

tt = class(Im);
%Im = cast(double(Im)/max(max(double(Im)))*gain,tt);
Im = cast(double(Im)*gain,tt);


[Imt,Imp]= remove_part(Im,part_radius);
% [Imt,Imp]= remove_part(Im,10);