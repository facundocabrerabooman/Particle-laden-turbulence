function CC = CenterFinding2D(folderin,folderout,fname,nstart,nframes)

%%
%N=7;
%R_range=[4 25];
%Thresh=40;
%sigma_low = 25;
%sigma_high = 2;

fnameMcin2 = [folderin filesep fname '.mcin2']

if isdeployed
    nframes=str2double(nframes);
    nstart=str2double(nstart);
end
%%



%% Find centers
folderout
if exist(folderout,'dir')==0
    mkdir(folderout);
end

fid = fopen([folderout filesep 'centers_' fname '.dat'],'w');
    fprintf(fid,'%d\n',nframes);

tic
for kframe=1:nframes
    
    [Im,Imref,params]=mCINREAD2(fnameMcin2,nstart+kframe-1,1);
%%
    Im=bsxfun(@minus,Im,Imref);
    Im=imclose(Im,strel('diamond',1));
    %Im=imerode(Im,strel('diamond',1));
    %Im(Im<=1)=0;
    %Imf = imfilter(Im,fspecial('gaussian',sigma_high*3,sigma_high),'symmetric');
    %Imf=imgaussfilt(Im,sigma_high,'Padding','symmetric');%-imgaussfilt(Im(:,:,1),sigma_low,'Padding','symmetric');
    ImM=imextendedmax(Im,1,4);
    %ImM=imregionalmax(Imf);
    PP=regionprops(ImM,'Centroid');
    cc=vertcat(PP.Centroid);
    x=cc(:,1);
    y=cc(:,2);
    CC(kframe).X=x';
    CC(kframe).Y=y';
    
 %%   
    
    fprintf(fid,'%d\n',kframe-1);
    fprintf(fid,'%d\n',numel(x));
    for kpart=1:numel(x);
        I=Im(round(y(kpart)),round(x(kpart)));
        Ncrop=3;
        Imc=imcrop(Im,[max(0,round(x(kpart)-Ncrop)) max(0,round(y(kpart)-Ncrop)) 2*Ncrop+1 2*Ncrop+1]);
        A=numel(find(Imc>I/3));
        CC(kframe).I(kpart)=I;
        CC(kframe).A(kpart)=A;
        fprintf(fid,'%6.3f %6.3f %d %d %6.3f %6.3f\n',x(kpart),y(kpart),I,I,A,A);
    end
    
%   find centers by 2 1D gaussian fits around maxima (x,y)

%     out=gauss2D_cntrd(Im,[round(x) round(y)],8);
%     CC(k).X = out(:,1);  % x center
%     CC(k).Y = out(:,2);  % y center
%     CC(k).Gax = out(:,3); % gaussian amplitude in x direction
%     CC(k).Gay = out(:,4); % gaussian amplitude in y direction
%     CC(k).Gsx = out(:,5); % gaussian width in x direction
%     CC(k).Gsy = out(:,6); % gaussian width in y direction
    
    
    
%     find centers by segmeting the image around maxima (x,y)
%
%     BW=imsegfmm(Im(:,:,1),round(x),round(y),.05); 
%     PP=regionprops('table',BW,'Centroid','FilledArea','Eccentricity');
%     
%     CC(k).X=PP.Centroid(:,1);
%     CC(k).Y=PP.Centroid(:,2);
%     CC(k).A=PP.FilledArea;
%     CC(k).E=PP.Eccentricity;
    
    
%   find centers by imfindcircles
%
%   [Ctmp,Rtmp] = imfindcircles(Imf,R_range,'EdgeThreshold',0.05,'Sensitivity',0.85);  
%     C(k).X=Ctmp(:,1);
%     C(k).Y=Ctmp(:,2);
%     C(k).R=Rtmp;
%     C(k).T=k;
end
toc
    fclose(fid);

% %% Track
% T=reshape(ones(1,N)'*(1:numel(C)),N*numel(C),[]);
% X=reshape([C.X],N*numel(C),[]);
% Y=reshape([C.Y],N*numel(C),[]);
% [tracks_lmin,tracks,t_start_stop]=track2d([X Y T],50,500);
% 
% %% convert tracks structure to vtracks
% Ntracks=max(tracks(:,4));
% for k=1:Ntracks
%     I=find(tracks(:,4)==k);
% 	vtracks(k).X=tracks(I,1);
%     vtracks(k).Y=tracks(I,2);
%     vtracks(k).T=tracks(I,3);
% end
