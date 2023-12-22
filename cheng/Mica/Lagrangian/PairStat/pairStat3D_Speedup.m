function [pair,trackpair]=pairStat3D_Speedup(part)
%part=pairStat(part)
%
% For a given part structure, calculates relative separation,
% velocity and acceleration of particles frame by frame.
% If you have a track strucuture, use track2part first.



% part=track2part(vtracks);
folder = pwd;
% disp(pwd)
%%
mkdir([folder '\temp_pair'])
parfor k=1:numel(part)
    save_temp_pair(folder,part,k)
end

%% sum up all files
disp('sum up ... temp_pair')
imax = 1e4;
pair(1:imax) = struct('dX',nan,'dY',nan,'dZ',nan,'dR2',nan,'dVx',nan,'dVy',nan,'dVz',nan,'dV2',nan,'dAx',nan,'dAy',nan,'dAz',nan,'dA2',nan,'dNtrack',nan);
for k=1:numel(part)
    load([folder '\temp_pair\temp_pair_' num2str(k) '.mat']);
    pair(k) = tp;
end
pair = pair(arrayfun(@(X)(~isnan(X.dX(1))),pair));
disp('deleting ... temp_pair')
rmdir([folder '\temp_pair'],'s')

% %% part2track for pairs
% dNtrack=[part.dNtrack];
% dNtrackU=unique(dNtrack);
% dX=[part.dX];
% dY=[part.dY];
% dR2=[part.dR2];
% dVx=[part.dVx];
% dVy=[part.dVy];
% dV2=[part.dV2];
% %dAx=[part.dAx];
% %dAy=[part.dAy];
% %dA2=[part.dA2];
% 
% for k=1:numel(dNtrackU)     
%    II=find(dNtrack==dNtrackU(k));  
%    trackpair(k).dX=dX(II);
%    trackpair(k).dY=dY(II);
%    trackpair(k).dR2=dR2(II);
%    trackpair(k).dVx=dVx(II);
%    trackpair(k).dVy=dVy(II);
%    trackpair(k).dV2=dV2(II);
%    trackpair(k).dNtrack=dNtrackU(k);
% end