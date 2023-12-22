clear
clc
close all

% input
fin = '';
fout = 'F:\matlaboutput\1.2g_1004_00V_3kHz_12dB';
%%
fout = [fout '\particle'];
if ~exist(fout, 'dir')
   mkdir(fout)
   disp(['make dir:' fout])
end

disp('read particle PTV')

dfolders = FunSubfolder(fin);

f = waitbar(0,'Please wait...');
for nexp = 1:size(dfolders,1)
    waitbar(nexp/size(dfolders,1), f, 'Please wait...');

    subfolders = [fin '\' dfolders(nexp).name];
    fname = [subfolders '\ExportToTecplot.dat'];
    [tracks,VARlist] = tec2mat(fname,'debug');

    tracklen(nexp) = size(tracks,2);

    Nframemax = 1e6;
    NTrackmax = 1e9;
    
    for in = 1:tracklen(nexp)
        numtracklen = sum(tracklen(1:nexp-1));
        NumPerFrame = size(tracks(in).data,1);
        if NumPerFrame>1
            warning([num2str(NumPerFrame) ' particles have been detected, the first one is saved '])
            disp(['Nexp = ' num2str(nexp) '; Nframe = ' num2str(str2double(tracks(in).T(end-3:end))+1)])
        end
        part.T(in+numtracklen,:) = (str2double(tracks(in).T(10:end))+Nframemax*(nexp-1))*ones(1,1);
        part.X(in+numtracklen,:) = tracks(in).data(1,1);
        part.Y(in+numtracklen,:) = tracks(in).data(1,2);
        part.Z(in+numtracklen,:) = tracks(in).data(1,3);
        part.Vx(in+numtracklen,:) = tracks(in).data(1,5);
        part.Vy(in+numtracklen,:) = tracks(in).data(1,6);
        part.Vz(in+numtracklen,:) = tracks(in).data(1,7);
        part.V(in+numtracklen,:) = tracks(in).data(1,8);
        part.Ntrack(in+numtracklen,:) = tracks(in).data(1,9)+NTrackmax*(nexp-1);
        part.Ax(in+numtracklen,:) = tracks(in).data(1,10);
        part.Ay(in+numtracklen,:) = tracks(in).data(1,11);
        part.Az(in+numtracklen,:) = tracks(in).data(1,12);
        part.A(in+numtracklen,:) = tracks(in).data(1,13);

    end
end


save([fout '\particle_PTV.mat'],'part')

waitbar(1, f, 'Done!');
pause(0.5)
close(f)