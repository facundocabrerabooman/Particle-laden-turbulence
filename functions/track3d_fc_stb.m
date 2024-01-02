function [traj,tracks]=track3d_fc_stb(data,folderin,folderout,fname,maxdist,lmin,flag_pred,npriormax,porder,flag_conf, numFrames, Fs)

% May 2022 FC
%
% Track particles on 3D images
% ____________________________________________________________________________
% INPUTS
% folderout : path to folder where track results will be saved
% fname     : name of movie files without extension (i.e. without .cam0.mcin2)
% maxdist   : maximum travelled distance between two successive frames
% lmin      : minimum length of a trajectory (number of frames)
% flag_pred : 1 for predictive tracking, 0 otherwise
% npriormax : maximum number of prior frames used for predictive tracking
% porder    : polynomial fitting order for predictive tracking
% flag_conf : 1 for conflict solving, 0 otherwise
%
% OUTPUTS
% traj(kt).ntraj  : trajectory index
% traj(kt).length : trajectory length
% traj(kt).frames : trajectory frames
% traj(kt).x      : x-position
% traj(kt).y      : y-position
% traj(kt).z      : z-position
% traj(kt).nmatch : element indices in tracks
% tracks          : trajectory raw data
% ____________________________________________________________________________

disp('Loading matches...');
tic

if isempty(data)
    load([folderin filesep fname '.mat'])
end
% if numFrames < size(d,1)
%     data=d(1:numFrames,:);
% else 
%     data=d;
% end
% clear d


toc 
disp('Tracking...');


%raw data
tracks=zeros(size(data,1),7);
tracks(:,1:4)=data(:,1:4); %frame number and position
tracks(:,5)=0; %trajectory index
tracks(:,6)=1; %availability of particles (0 unavailable, 1 available, 2 conflict)
tracks(:,7)=(1:length(data(:,1)))'; %element index

%number of particles in first frame = number of active trajectories
ind_act=find(tracks(:,1)==1);               % fcb changed this from ==1 to ==2
tracks(ind_act,5)=ind_act;
tracks(ind_act,6)=0;
ntraj=max(ind_act);

%length of trajectories (easy elimination of short trajectories later)
lmax=1e5;
ltraj=zeros(lmax,2);
ltraj(1:ntraj,1)=(1:ntraj)'; %trajectory index
ltraj(1:ntraj,2)=1; %each trajectory starts with one particle

for kf=2:max(tracks(:,1)) %for each frame
    
    ind_new=find(tracks(:,1)==kf);
    min_new=min(ind_new);
    max_new=max(ind_new);
    c=1;
    dispmax=1e5;
    dispoind=zeros(dispmax,1);
    dispoact=zeros(dispmax,1);
    disposqrd=zeros(dispmax,1);
    
    %positions of particles in frame kf
    newx=tracks(ind_new,2);
    newy=tracks(ind_new,3);
    newz=tracks(ind_new,4);

    for kp=1:length(ind_act) %for each particle in frame kf-1
        
        %particle kp in frame kf-1
        actx=tracks(ind_act(kp),2); %x-position
        acty=tracks(ind_act(kp),3); %y-position
        actz=tracks(ind_act(kp),4); %z-position 
        actnum=tracks(ind_act(kp),5); %trajectory index
        
        %predictive tracking
        if flag_pred==1
            %indices of up to last npriormax points in active trajectory
            ind_pred=find(tracks(:,5)==actnum,npriormax,'last');
            %length of used active trajectory part
            nprior=length(ind_pred);
            if nprior>1
                porder_eff=min(porder,nprior-1); %degree < number of points
                %fit trajectory
                fitx=polyfit((1:nprior)',tracks(ind_pred,2),porder_eff);
                fity=polyfit((1:nprior)',tracks(ind_pred,3),porder_eff);
                fitz=polyfit((1:nprior)',tracks(ind_pred,4),porder_eff);
                %extrapolate position in next frame
                actx=polyval(fitx,nprior+1); %x-position
                acty=polyval(fity,nprior+1); %y-position
                actz=polyval(fitz,nprior+1); %z-position
            end
        end
        
        %calculation of squared distance: nearest neighbour
        sqrdist=(actx-newx).^2+(acty-newy).^2+(actz-newz).^2;
        
        %minimize the distance
        [sqrdmin,ind_min]=min(sqrdist);
        if sqrdmin<maxdist^2
            dispo=tracks(ind_new(ind_min),6);
            
            %available particle
            if dispo==1 
                %particle becomes unavailable
                tracks(ind_new(ind_min),6)=0;
                %new trajectory index
                tracks(ind_new(ind_min),5)=actnum;
                %lengthen the trajectory
                ltraj(actnum,2)=ltraj(actnum,2)+1;
                %useful arrays for next steps
                dispoind(c)=ind_new(ind_min); %element index
                dispoact(c)=actnum; %trajectory index
                disposqrd(c)=sqrdmin; %minimum squared distance
                c=c+1;
                
            %unavailable particle = conflict
            elseif dispo==0
                i=find(dispoind==ind_new(ind_min));
                actnump=dispoact(i);
                
                %non solving
                if flag_conf==0 || sqrdmin==disposqrd(i)
                    %particle could go to two different trajectories
                    tracks(ind_new(ind_min),6)=2;
                    %trajectory index set to zero 
                    tracks(ind_new(ind_min),5)=0;
                    %shorten the prior trajectory
                    ltraj(actnump,2)=ltraj(actnump,2)-1;
                end
                
                %solving
                if flag_conf==1 && sqrdmin<disposqrd(i) %new track is a better match
                    %new trajectory index
                    tracks(ind_new(ind_min),5)=actnum;
                    %lengthen the trajectory
                    ltraj(actnum,2)=ltraj(actnum,2)+1;
                    %shorten the prior trajectory
                    ltraj(actnump,2)=ltraj(actnump,2)-1;
                    %array updating
                    dispoact(i)=actnum;
                    disposqrd(i)=sqrdmin;
                end

            end
            
        end
        
    end
    
    %tracked particules: found once or not found
    tracks_temp=tracks(min_new:max_new,:);
    %found once
    ind=tracks_temp(:,6)==0;
    ind_act0=tracks_temp(ind,7);
    %not found
    ind=tracks_temp(:,6)==1;
    ind_act1=tracks_temp(ind,7);
    if ~isempty(ind_act1) %new trajectory index for particles not found
        for k=1:length(ind_act1)
            ntraj=ntraj+1;
            tracks(ind_act1(k),5)=ntraj;
            ltraj(ntraj,1)=ntraj;
            ltraj(ntraj,2)=1;
        end
%     elseif  isempty(ntraj)% fcb may9 23
%         ntraj=ntraj;
%             tracks(ind_act1(k),5)=[];
%             ltraj(ntraj,1)=[];
%             ltraj(ntraj,2)=[];
    end
    ind_act=[ind_act0;ind_act1];
    
    %display frame number every 100 frames
    if rem(kf,100)==0 
        disp(['Trajectories found for ' num2str(kf) ' frames']);
    end
    
end

%write results
i0=tracks(:,5); %index of trajectories
i=find(ltraj(:,2)>=lmin); %index of trajectories >= lmin
traj=struct([]);
for kt=1:length(i) %for each trajectory
    I=find(i0==i(kt));
    lt=length(I);
    traj(kt).ntraj=kt; %trajectory index
    traj(kt).length=lt; %trajectory length
    traj(kt).Ntrack=tracks(I,1); %trajectory frames
    traj(kt).t=tracks(I,1)./Fs; %trajectory times
    traj(kt).x=tracks(I,2); %x-position
    traj(kt).y=tracks(I,3); %y-position
    traj(kt).z=tracks(I,4); %z-position
    traj(kt).nmatch=tracks(I,7); %element indices
end

disp([num2str(kt) ' trajectories longer than ' num2str(lmin) ' frames (from ' num2str(length(vertcat(traj.x))) ' matches)']);
disp(['Saving to .mat file in ' folderout]);
save([folderout filesep 'trajs_' fname '.mat'],'traj','tracks','-v7.3');
save([folderout filesep 'param_trajs_' fname '.mat'],'maxdist','lmin','flag_pred','npriormax','porder','flag_conf')
disp('Tracking complete!');
elapsed=toc;
disp(['Elapsed time: ',num2str(elapsed,'%.2f'),'s']);

 end