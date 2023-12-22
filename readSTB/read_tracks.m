clear
clc
close all
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'))

% input
fin = '\\nasvierstokes\home\Cheng\workspace_lavision_00V\1.0g_0929_00V_3kHz_12dB\preproc_particle\AddGeometricMask\PTV\Repair particle tracks\ExportToTecplot';
fout = '\\nasvierstokes\home\Cheng\output_matlab_00V\1.0g_0929_00V_3kHz_12dB';

FlagInertial = 'on';

%% save to seperate files
if strcmp(FlagInertial,'on')
    fout = [fout '\particle'];
elseif strcmp(FlagInertial,'off')
    fout = [fout '\tracers'];
end

if ~exist(fout, 'dir')
   mkdir(fout)
   disp(['make dir:' fout(end-6:end)])
end

disp('read Davis Output')

dfolders = FunSubfolder(fin);

%%
if size(dfolders,1) ~= 0
    if ~exist([fout '\Seperate_Files_STB'],'dir')
        mkdir([fout '\Seperate_Files_STB'])
        disp('make dir: Seperate_Files_STB')
    end

    parfor nexp = 1:size(dfolders,1)
        subfolders = [fin '\' dfolders(nexp).name];
        fname = [subfolders '\ExportToTecplot.dat']; % check the file name, could be: ShakeTheBox tr.dat
        
        stb2mat(fname,fout,FlagInertial,nexp)
    end
else
    fname = [fin '\ExportToTecplot.dat'];
    stb2mat(fname,fout,FlagInertial)
end
disp('seperate files saved!')

%% sum up
if size(dfolders,1) ~= 0
    temp_part = [];
    f = waitbar(0,'Please wait...');
        for i = 1:size(dfolders,1)
            waitbar(i/size(dfolders,1), f, 'Please wait...');
            fnamet = [fout '\Seperate_Files_STB\STB_' num2str(i) '.mat'];
            if exist(fnamet) 
                t = load(fnamet);
                temp_part = [temp_part;t.part];
            end
        end
        fld1 = fieldnames(temp_part);     
        for ii = 1:length(fld1)
            part.(fld1{ii}) = vertcat(temp_part.(fld1{ii}));
        end
    waitbar(1, f, 'Saving...');
    save([fout '\STB.mat'],'part', '-v7.3')
    waitbar(1, f, 'Done!');
    pause(0.5)
    close(f)
end


%% function definition
function stb2mat(fname,fout,FlagInertial,n) 
    if nargin<4
        n=1; 
    end
    [tracks,~] = tec2mat(fname,'debug');

    tracklen = size(tracks,2);

    Nframemax = 1e6;
    NTrackmax = 1e9;

    part = [];

    for nf = 1:tracklen
        if strcmp(FlagInertial,'on')
            NumPerFrame(nf) = size(tracks(nf).data,1);
            if NumPerFrame(nf)>1
                warning([num2str(NumPerFrame(nf)) ' particles have been detected, the first one is saved '])
                disp(['Nexp = ' num2str(n) '; Nframe = ' num2str(str2double(tracks(nf).T(end-3:end))+1)])
                NumPerFrame(nf) = 1;
            end
            ind = (sum(NumPerFrame(1:nf-1))+1):1:sum(NumPerFrame);
            part.T(ind,:) = (str2double(tracks(nf).T(10:end))+Nframemax*(n-1))*ones(NumPerFrame(nf),1);
            part.Ntrack(ind,:) = tracks(nf).data(1,9)+NTrackmax*(n-1);
            part.X(ind,:) = tracks(nf).data(1,1);
            part.Y(ind,:) = tracks(nf).data(1,2);
            part.Z(ind,:) = tracks(nf).data(1,3);
            part.Vx(ind,:) = tracks(nf).data(1,5);
            part.Vy(ind,:) = tracks(nf).data(1,6);
            part.Vz(ind,:) = tracks(nf).data(1,7);
            part.V(ind,:) = tracks(nf).data(1,8);
            part.Ax(ind,:) = tracks(nf).data(1,10);
            part.Ay(ind,:) = tracks(nf).data(1,11);
            part.Az(ind,:) = tracks(nf).data(1,12);
            part.A(ind,:) = tracks(nf).data(1,13);
        elseif strcmp(FlagInertial,'off')
            NumPerFrame(nf) = size(tracks(nf).data,1);
            ind = (sum(NumPerFrame(1:nf-1))+1):1:sum(NumPerFrame);
            part.T(ind,:) = (str2double(tracks(nf).T(10:end))+Nframemax*(n-1))*ones(NumPerFrame(nf),1);
            part.Ntrack(ind,:) = tracks(nf).data(:,9)+NTrackmax*(n-1);
            part.X(ind,:) = tracks(nf).data(:,1);
            part.Y(ind,:) = tracks(nf).data(:,2);
            part.Z(ind,:) = tracks(nf).data(:,3);
        end
    end
    if ~isempty(part) % in case particle is not found
        % save for each folder to speed up the processing
        if nargin>3
            save([fout '\Seperate_Files_STB\STB_' num2str(n) '.mat'],'part')
        else
            save([fout '\STB.mat'],'part')
        end
    end
end