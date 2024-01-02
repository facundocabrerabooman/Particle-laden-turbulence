function [sliceData,subX, subY, subZ] = extractSlice_Cheng(grids,fields,center,normVector,radius,nSpacing)


% if centerX,centerY,centerZ = 0, subX,subY,subZ is the same as
% sliceData.XData,YData,ZData.



if nargin < 4
   display('requires at least 7 parameters');
   return;
end

if nargin < 5
    % sets the size for output slice radius*2+1.
    radius = 10;
end

pt = [center.x,center.y,center.z];
vec = [normVector.x,normVector.y,normVector.z];

%initialize needed parameters
%size of volume.
volSz=size(fields); 
%a very small value.
epsilon = 1e-12; 

%assume the slice is initially at [0,0,0] with a vector [0,0,1] and a silceSize.
% hsp = surf(linspace(-radius,radius,2*radius+1),linspace(-radius,radius,2*radius+1),zeros(2*radius+1));
hsp = surf(linspace(-radius,radius,nSpacing),linspace(-radius,radius,nSpacing),zeros(nSpacing));
hspInitialVector = [0,0,1];

%normalize vectors;
hspVec = hspInitialVector/norm(hspInitialVector);
hspVec(hspVec ==0) = epsilon;
vec = vec/norm(vec);
vec(vec == 0)=epsilon;

%this does not rotate the surface, but initializes the subscript z in hsp.
rotate(hsp,[0,0,1],0);

%get the coordinates
xdO = get(hsp,'XData');
ydO = get(hsp,'YData');
zdO = get(hsp,'ZData');

hspVecXvec = cross(hspVec, vec)/norm(cross(hspVec, vec));
acosineVal = acos(dot(hspVec, vec));

%help prevents errors (i.e. if vec is same as hspVec),
hspVecXvec(isnan(hspVecXvec)) = epsilon;
acosineVal(isnan(acosineVal)) = epsilon;

%rotate to the requested orientation
rotate(hsp,hspVecXvec(:)',180*acosineVal/pi);

%get the coordinates
xd = get(hsp,'XData');
yd = get(hsp,'YData');
zd = get(hsp,'ZData');


% figure;
% plot3(xdO,ydO,zdO,'b+'); hold on;
% plot3(xd,yd,zd,'m.'); hold off

%translate;
subX = xd + pt(1);
subY = yd + pt(2);
subZ = zd + pt(3);


%%
sdata = slice(grids.XX,grids.YY,grids.ZZ,fields,xd,yd,zd,'linear');
sliceData.xData = sdata.XData;
sliceData.yData = sdata.YData;
sliceData.zData = sdata.ZData;
sliceData.cData = sdata.CData;
close
% sliceData.xData = xd; sliceData.yData = yd; sliceData.zData = zd;
