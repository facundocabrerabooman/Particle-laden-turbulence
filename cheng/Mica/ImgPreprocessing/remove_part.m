function [Imt,Imp] = remove_part(Im,PartRadius,thres_part,iffilter)

st = strel('disk',PartRadius);
stdil = strel('disk',round(PartRadius/3));

Imd = imdilate(imbinarize(imopen(Im,st)),stdil);

%% refine particle image
Rp=sqrt(numel(find(Imd~=0)))/pi;

Rmin = max(0.6*Rp,5.5);
Rmax = max(1.6*Rp,15);
[C, R] = imfindcircles(Imd, cast([Rmin Rmax ],class(Im)));
if size(C,1) == 1
    roi=images.roi.Circle('Center',C,'Radius',R);
    mask = createMask(roi,size(Imd,1),size(Imd,2));
else
    mask = Imd;
end
Imp = immultiply(Im,cast(mask,class(Im)));
Imt = immultiply(Im,cast(abs(1-mask),class(Im)));

Imp(Imp<thres_part)=0;

Imp = imgaussfilt(Imp,1);
% Imt = imgaussfilt(Imt,1);
if iffilter==1
%     Imp = imgaussfilt(Imp,1);
    Imt = imgaussfilt(Imt,1);
end














