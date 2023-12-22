clear
clc
close all

% input
fin = 'F:\LavisionWorkspace\0811_PinholeCalib_1g_12V_4kfps\tracers_notNormalized\ShakeTheBox tr\ExportToTecplot';

fout = 'F:\MaltlabOutput\0811_1g_12V_4kfps';
%% save to seperate files
fout = [fout '\tracers'];
if ~exist(fout, 'dir')
   mkdir(fout)
   disp(['make dir:' fout])
end

disp('read tracers STB')

dfolders = FunSubfolder(fin);
if size(dfolders,1) ~= 0
    parfor nexp = 1:size(dfolders,1)
        subfolders = [fin '\' dfolders(nexp).name];
        fname = [subfolders '\ExportToTecplot.dat']; % check the file name, could be: ShakeTheBox tr.dat
        
        save2tracerSTB(fname,fout,nexp)
    end
else
    fname = [fin '\ShakeTheBox tr.dat'];
    save2tracerSTB(fname,fout)
end
disp('seperate files saved!')

%% sum up
if size(dfolders,1) ~= 0
    tracer = [];
    f = waitbar(0,'Please wait...');
        for i = 1:size(dfolders,1)
            waitbar(i/size(dfolders,1), f, 'Please wait...');
            t = load([fout '\tracer_STB_' num2str(i) '.mat']);
            tracer = [tracer,t.tracer];
        end
    waitbar(1, f, 'Saving...');
    save([fout '\tracer_STB.mat'],'tracer', '-v7.3')
    waitbar(1, f, 'Done!');
    pause(0.5)
    close(f)
end


%% function definition
function save2tracerSTB(fname,fout,n) 
    if nargin<3
        n=1; 
    end
    [tracks,~] = tec2mat(fname,'debug');

    tracklen = size(tracks,2);

    Nframemax = 1e6;
    NTrackmax = 1e9;
    
    fmax = 1e5; % max frames 
    tracer(1:fmax) = struct('T',nan,'X',nan,'Y',nan,'Z',nan,'Ntrack',nan);

    for nf = 1:tracklen
%         NumPerFrame = size(tracks(nf).data,1);
        %tracer(nf).Nexp = n;
        tracer(nf).T = (str2double(tracks(nf).T(10:end))+Nframemax*(n-1))*ones(size(tracks(nf).data,1),1);
        tracer(nf).X = tracks(nf).data(:,1);
        tracer(nf).Y = tracks(nf).data(:,2);
        tracer(nf).Z = tracks(nf).data(:,3);
        %tracer(nf).vx = tracks(nf).data(:,5);
        %tracer(nf).vy = tracks(nf).data(:,6);
        %tracer(nf).vz = tracks(nf).data(:,7);
        %tracer(nf).V = tracks(nf).data(:,8);
        tracer(nf).Ntrack = tracks(nf).data(:,9)+NTrackmax*(n-1);
        %tracer(nf).ax = tracks(nf).data(:,10);
        %tracer(nf).ay = tracks(nf).data(:,11);
        %tracer(nf).az = tracks(nf).data(:,12);
        %tracer(nf).A = tracks(nf).data(:,13);
    end
    tracer = tracer(arrayfun(@(X)(~isnan(X.T(1))),tracer));
    % save for each folder to speed up the processing
    if nargin>2   
        save([fout '\Seperate_Files_STB\tracer_STB_' num2str(n) '.mat'],'tracer')
    else
        save([fout '\tracer_STB.mat'],'tracer')
    end
end