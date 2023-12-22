function RetainStruct = retainNeighborTracks(part,tracer,outStructFields,idx_front_back)

fieldname = fieldnames(tracer);
RetainData = struct();

%% copy the neighbor tracks
% coordinates need to be convert to the particle framework
% velocity and acceleration need to minus the particle's V and A
% rotation matrix needs to be applied to align the particle's moving 
% direction cross different frames and experiments. The rotation matrix is
% built by quaternion. We align the particle motion to verticle

% vector of the gravity, be careful the switch of data in X and Y axis
% notes: in part and tracers (trajectory structure), X and Z are horizontal
% component, Y is the vertical direction (gravity)
vectorg = [0,-1,0];

for i = 1:numel(idx_front_back)
    
    %% index for T, Ntrack
    idxp = find(part.T == tracer.Tf(idx_front_back(i)));
    idx0a = find(tracer.T == part.T(idxp));
    idx0b = tracer.Ntrack(idx0a) == tracer.Ntrackf(idx_front_back(i));
    idxt = idx0a(idx0b);
    
    RetainData.T(i,1) = tracer.T(idxt);
    RetainData.Ntrack(i,1) = tracer.Ntrack(idxt);

    %% get the particle's velocity vector of the current instant time,
    % and prepare the rotation matrix
    idxp = part.Tf == tracer.Tf(idx_front_back(i));
    Vp = [part.Vx(idxp); part.Vy(idxp); part.Vz(idxp)];
    R = rotationMatrix(Vp, vectorg);
    
    %% index for X,Y,Z
    idxp = find(part.T == tracer.Tf(idx_front_back(i)));
    idx0a = find(tracer.T == part.T(idxp));
    idx0b = tracer.Ntrack(idx0a) == tracer.Ntrackf(idx_front_back(i));
    idxt = idx0a(idx0b);
    
    coordParticleFramework(1,1) = tracer.X(idxt)-part.X(idxp);
    coordParticleFramework(2,1) = tracer.Y(idxt)-part.Y(idxp);
    coordParticleFramework(3,1) = tracer.Z(idxt)-part.Z(idxp);

    coordRotated = R*coordParticleFramework;
    RetainData.X(i,1) = coordRotated(1,1);
    RetainData.Y(i,1) = coordRotated(2,1);
    RetainData.Z(i,1) = coordRotated(3,1);

    %% index for Xf,Yf,Zf,Vx,Vy,Vz,Ax,Ay,Az
    idxp = part.Tf == tracer.Tf(idx_front_back(i));
    idxt = idx_front_back(i);

    coordfParticleFramework(1,1) = tracer.Xf(idxt)-part.Xf(idxp);
    coordfParticleFramework(2,1) = tracer.Yf(idxt)-part.Yf(idxp);
    coordfParticleFramework(3,1) = tracer.Zf(idxt)-part.Zf(idxp);
    velofParticleFramework(1,1) = tracer.Vx(idxt)-part.Vx(idxp);
    velofParticleFramework(2,1) = tracer.Vy(idxt)-part.Vy(idxp);
    velofParticleFramework(3,1) = tracer.Vz(idxt)-part.Vz(idxp);
    accefParticleFramework(1,1) = tracer.Ax(idxt)-part.Ax(idxp);
    accefParticleFramework(2,1) = tracer.Ay(idxt)-part.Ay(idxp);
    accefParticleFramework(3,1) = tracer.Az(idxt)-part.Az(idxp);

    coordfRotated = R*coordfParticleFramework;
    velofRotated = R*velofParticleFramework;
    accefRotated = R*accefParticleFramework;

    RetainData.Xf(i,1) = coordfRotated(1,1);
    RetainData.Yf(i,1) = coordfRotated(2,1);
    RetainData.Zf(i,1) = coordfRotated(3,1);
    RetainData.Vx(i,1) = velofRotated(1,1);
    RetainData.Vy(i,1) = velofRotated(2,1);
    RetainData.Vz(i,1) = velofRotated(3,1);
    RetainData.Ax(i,1) = accefRotated(1,1);
    RetainData.Ay(i,1) = accefRotated(2,1);
    RetainData.Az(i,1) = accefRotated(3,1);
       
    %% index for Tf, Ntrackf
    idxt = idx_front_back(i);
    RetainData.Tf(i,1) = tracer.Tf(idxt);
    RetainData.Ntrackf(i,1) = tracer.Ntrackf(idxt);
   
%     disp([num2str(idxp) '   ' fieldname{j} '   ' 'tracer' num2str(tracer.(fieldname{j})(idxt)) ' '   'part' num2str(part.(fieldname{j})(idxp))])
end

%% convert to tracjectory structure
ntrack = unique(RetainData.Ntrack);
for i = 1:numel(ntrack)
    idxp = RetainData.Ntrack == ntrack(i);
    for j = 1:numel(outStructFields)
        RetainStruct(i,1).(outStructFields{j}) = RetainData.(outStructFields{j})(idxp);
    end
end