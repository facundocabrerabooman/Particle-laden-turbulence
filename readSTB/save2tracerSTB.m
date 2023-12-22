function save2tracerSTB(fname,fout,n)
    
    [tracks,~] = tec2mat(fname,'debug');

    tracklen = size(tracks,2);

    Nframemax = 1e6;
    NTrackmax = 1e9;

    for nf = 1:tracklen
        NumPerFrame = size(tracks(nf).data,1);
        %tracer(nf).Nexp = n;
        tracer(nf).T = (str2double(tracks(nf).T(end-3:end))+Nframemax*(n-1))*ones(size(tracks(nf).data,1),1);
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
    
    % save for each folder to speed up the processing
    save([fout '\tracer_STB_' num2str(n) '.mat'],'tracer')