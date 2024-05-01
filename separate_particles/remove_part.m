function [Imt,Imp] = remove_part(Im,PartRadius)

st = strel('disk',PartRadius);
%stdil = strel('disk',round(PartRadius/3));

%Imd = imdilate(imbinarize(imopen(Im,st)),stdil);
Imd = imbinarize(imopen(Im,st));
%Imd=Im;
%% refine particle image
%Rp=sqrt(numel(find(Imd~=0)))/pi;

Rmin = max(0.4*PartRadius,5);
Rmax = max(1.8*PartRadius,15);

Rmin = max(0.8*PartRadius,5);
Rmax = max(1.2*PartRadius,15);

[C, R] = imfindcircles(Imd, cast([Rmin Rmax],class(Im)));
%if size(C,1) == 1
%    roi=images.roi.Circle('Center',C,'Radius',R);
%    mask = createMask(roi,size(Imd,1),size(Imd,2));
%    Imp = immultiply(Im,cast(mask,class(Im)));
%Imt = immultiply(Im,cast(abs(1-mask),class(Im)));  
%else
    %Imt = Imd;
    %Imp = zeros(size(Imd,1),size(Imd,2));
%end
if size(C,1) == 1
    roi=images.roi.Circle('Center',C,'Radius',R);
    mask = createMask(roi,size(Imd,1),size(Imd,2));
else
    mask = Imd;
end
Imp = immultiply(Im,cast(mask,class(Im)));
Imt = immultiply(Im,cast(abs(1-mask),class(Im)));


Imp = imgaussfilt(Imp,1);
Imt = imgaussfilt(Imt,1);













