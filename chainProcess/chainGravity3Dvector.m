clc;clear;close all;

% input
fin = 'G:\LavisionWorkspace\calib_1016\chain_preproc\9x9 smoothing\PTV\ExportToTecplot';
fout = 'G:\matlaboutput\calib_1016';

tic
%%
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));
if ~exist(fout, 'dir')
   mkdir(fout)
   disp(['make dir:' fout])
end

disp('read chain PTV')


f = waitbar(0,'Please wait...');

fname = [fin '\PTV.dat'];
[tracks,VARlist] = tec2mat(fname,'debug');

tracklen = size(tracks,2);

for nf = 1:tracklen
    part(nf).T = str2double(tracks(nf).T(10:end))+zeros(size(tracks(nf).data,1),1);
    part(nf).X = tracks(nf).data(:,1);
    part(nf).Y = tracks(nf).data(:,2);
    part(nf).Z = tracks(nf).data(:,3);
    part(nf).Ntrack = tracks(nf).data(:,9);
end

track = part2track(part);

for i = 1:size(track,2)
    xyz(1,i) = mean(track(i).X);
    xyz(2,i) = mean(track(i).Y);
    xyz(3,i) = mean(track(i).Z);
end

%%
% https://fr.mathworks.com/matlabcentral/answers/424591-3d-best-fit-line
% xzyl is 3x2 contains coordinates of 2 points of the lines
% xzyl(:,1) is the 3D coordinates of the first point
% xzyl(:,2) is the 3D coordinates of the second point
% The line equation (parametric form) is then
% xyz(t) = t*xzyl(:,1) + (1-t)*(xzyl(:,2);
% where t is a real parameter.

xyz0 = mean(xyz,2);
A = xyz-xyz0;
[U,S,~] = svd(A);
d = U(:,1);
t = d'*A;
t1 = min(t);
t2 = max(t);
xzyl = xyz0 + [t1,t2].*d; % size 3x2
% Check
x = xyz(1,:);
y = xyz(2,:);
z = xyz(3,:);
xl = xzyl(1,:);
yl = xzyl(2,:);
zl = xzyl(3,:);
close all

plot3(x,y,z,'o');
hold on
plot3(xl,yl,zl,'r');
axis equal


%%
save([fout '\chain_PTV.mat'],'xyz','xzyl')

waitbar(1, f, 'Done!');
pause(0.5)
close(f)
toc