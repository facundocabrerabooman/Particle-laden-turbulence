function RetainStruct = retainNeighborTracks_old(part,tracer,outStructFields,idx_front_back)

fieldname = fieldnames(tracer);
RetainData = struct();

%% copy the neighbor tracks
% coordinates need to be convert to the particle framework
% velocity and acceleration need to minus the particle's V and A
% rotation matrix needs to be applied to align the particle's moving 
% direction cross different frames and experiments. The rotation matrix is
% built by quaternion. We align the particle motion to verticle

% vector of the gravity, be careful the swithch of data in X and Y axis
vectorg = [-1,0,0];

for i = 1:numel(idx_front_back)
    for j = 1:numel(fieldname)
        if j<=2 % T, Ntrack
            idxp = find(part.T == tracer.Tf(idx_front_back(i)));
            idx0a = find(tracer.T == part.T(idxp));
            idx0b = tracer.Ntrack(idx0a) == tracer.Ntrackf(idx_front_back(i));
            idxt = idx0a(idx0b);
            RetainData.(fieldname{j})(i,1) = tracer.(fieldname{j})(idxt);
        elseif j>=3 && j<=5 % X,Y,Z
            idxp = find(part.T == tracer.Tf(idx_front_back(i)));
            idx0a = find(tracer.T == part.T(idxp));
            idx0b = tracer.Ntrack(idx0a) == tracer.Ntrackf(idx_front_back(i));
            idxt = idx0a(idx0b);
            RetainData.(fieldname{j})(i,1) = tracer.(fieldname{j})(idxt)-part.(fieldname{j})(idxp);
        elseif j>=6 && j<=14 % Xf,Yf,Zf, vx,vy,vz,ax,ay,az
            idxp = part.Tf == tracer.Tf(idx_front_back(i));
            idxt = idx_front_back(i);
            RetainData.(fieldname{j})(i,1) = tracer.(fieldname{j})(idxt)-part.(fieldname{j})(idxp);
%         elseif j>=9 && j<=14 % vx,vy,vz,ax,ay,az
%             idxp = find(part.Tf == tracer.Tf(idx_front_back(i)));
%             idxt = idx_front_back(i);
%             RetainData.(fieldname{j})(i,1) = tracer.(fieldname{j})(idxt);
        else % Tf, Ntrackf
%             idxp = find(part.Tf == tracer.Tf(idx_front_back(i)));
            idxt = idx_front_back(i);
            RetainData.(fieldname{j})(i,1) = tracer.(fieldname{j})(idxt);
        end
%         disp([num2str(idxp) '   ' fieldname{j} '   ' 'tracer' num2str(tracer.(fieldname{j})(idxt)) ' '   'part' num2str(part.(fieldname{j})(idxp))])
    end
end

%% convert to tracjectory structure
ntrack = unique(RetainData.Ntrack);
for i = 1:numel(ntrack)
    idxp = RetainData.Ntrack == ntrack(i);
    for j = 1:numel(outStructFields)
        RetainStruct(i,1).(outStructFields{j}) = RetainData.(outStructFields{j})(idxp);
    end
end