function centers=centerfinding(folderin,folderout,fname,nstart,nframes)

% centers=centerfinding(folderin,folderout,fname,nstart,nframes)
% 03/2019 - Thomas Basset
%
% Find centers of particles on 2D images
% ____________________________________________________________________________
% INPUTS
% folderin  : path to folder containing movies
% folderout : path to folder where center results will be saved
% fname     : name of movie files without extension (i.e. without .cam0.mcin2)
% nstart    : first frame to be loaded
% nframes   : number of frames to be loaded
%
% OUTPUTS
% centers(kcam).frame(kf).x     : x-position
% centers(kcam).frame(kf).y     : y-position
% centers(kcam).frame(kf).majax : major axis length
% centers(kcam).frame(kf).minax : minor axis length
% centers(kcam).frame(kf).A     : area
% centers(kcam).frame(kf).I     : intensity
% ____________________________________________________________________________

disp(['First frame read: ' num2str(nstart)]);
disp(['Number of frames read: ' num2str(nframes)]);
disp('Beginning center finding...');
tic

for kcam=0:2 %for each camera (use parfor to parallelize)
    
    %input file name
    movie=[folderin filesep fname '.cam' num2str(kcam) '.mcin2'];
    disp(['Input file: ' movie]);
    
    %read frames
    [im0,imref,~]=mCINREAD2(movie,nstart,nframes);
    
    %subtract black reference
    imref=repmat(imref,1,1,nframes);
    im=imsubtract(im0,imref); 
    
    %subtract mean reference
    all_im=mCINREAD2(movie);
    imrefmean=mean(all_im,3);
    imrefmean=repmat(imrefmean,1,1,nframes);
    im=imsubtract(double(im),imrefmean); 
    im=imcomplement(im);
    im=uint8(im);
    
    %detect particles
    frame=struct([]);
    for kf=1:nframes %for each frame (use parfor to parallelize)
        Im=im(:,:,kf);
        c=findcenters_tb(Im);
        frame(kf).x=c(:,1);
        frame(kf).y=c(:,2);
        frame(kf).majax=c(:,3);
        frame(kf).minax=c(:,4);
        frame(kf).A=c(:,5);
        frame(kf).I=c(:,6);
        %display frame number every 100 frames
        if rem(kf,100)==0 
            disp(['Centers found for ' num2str(kf) ' frames']);
        end
    end
    centers(kcam+1).frame=frame;

end

%output folder name
disp(['Saving to .dat files in ' folderout]);
%if directory doesn't exist, create one
if exist(folderout,'dir')==0
    mkdir(folderout);
end

%write results
for kcam=0:2 %for each camera (use parfor to parallelize)
    fid=fopen([folderout filesep 'centers_' fname '.cam' num2str(kcam) '.dat'],'w');
    fprintf(fid,'%d\n',nframes); %print number of frames
    for kf=1:nframes %for each frame
        fprintf(fid,'%d\n',kf); %print frame number
        npart=length(centers(kcam+1).frame(kf).x);
        fprintf(fid,'%d\n',npart); %print number of detections
        for kp=1:npart %for each particle
        fprintf(fid,'%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',...
            centers(kcam+1).frame(kf).x(kp),centers(kcam+1).frame(kf).y(kp),centers(kcam+1).frame(kf).majax(kp),...
            centers(kcam+1).frame(kf).minax(kp),centers(kcam+1).frame(kf).A(kp),centers(kcam+1).frame(kf).I(kp)); %print data
        end
    end
    fclose(fid);
end

disp('Center finding complete!');
elapsed=toc;
disp(['Elapsed time: ' num2str(elapsed,'%.2f') 's']);
    
end