function [MeanSlice,sliceData,varargout]  = meanSlice(grids,meanField,nSlice,Rmax,nSpacing)


initial_normal_vector = [0, 0, 1]; 
rotation_axis = 'X'; 
% rotation_angle = 3*pi/2; 
angle_rng = pi;


center.x = 0;
center.y = 0;
center.z = 0;
radius = Rmax;


if center.x ~= 0 || center.y ~= 0 || center.z ~= 0
    varargout = struct('subX',[],'subY',[],'subZ',[]);
end

for i = 1:nSlice
    rotation_angle = (i-1)/nSlice*angle_rng;
    normVector = rotate_plane(initial_normal_vector, rotation_axis, rotation_angle);
    
    [sliceData(i),subX, subY, subZ]= extractSlice_Cheng(grids,meanField,center,normVector,radius,nSpacing);

%     surf(subX,subY,subZ,nslice,'FaceColor','texturemap','EdgeColor','none');hold on
%     drawnow;

    if center.x ~= 0 || center.y ~= 0 || center.z ~= 0
        varargout(i).subX = subX;
        varargout(i).subY = subY;
        varargout(i).subZ = subZ;
    end

end


MeanSlice = ones(size(sliceData(1).xData))*NaN;
for i = 1:size(sliceData(1).xData,1)
    for j = 1:size(sliceData(1).xData,2)
        MeanSlice(i,j) = mean(arrayfun(@(x) vertcat(x.cData(i,j)),sliceData),'omitnan');
    end
end