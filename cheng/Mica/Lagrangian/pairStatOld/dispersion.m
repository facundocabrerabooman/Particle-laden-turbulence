function [X1,X2,Y1,Y2,Z1,Z2,d2,len,nt]=dispersion(minoverlap,file0,num_of_files,filename_save,varargin)

% 15/11/2004
% analyze the dispersion of pairs from the track files
%

%load tracks

%frame0[i]=first frame of track i
%nframes[i]=number of frame in track i
%maxlength[i]=length of longest track

%1st find the overlap tensor O
%
%   Oij = number of frames common to track i and j~=i
%       =  0 if i=j

if(nargin>4)
    filter_w=varargin{1};
    filter_l=varargin{2};
    ker=posfiltcoef(filter_w,filter_l);
    kerv=velfiltcoef(filter_w,filter_l);
else
    filter_w=0;
    filter_l=0;
end

minlen=minoverlap;
maxoverlap=2000;
nt=zeros(1,maxoverlap);
%minoverlap=200;
ii_tot=0;
for nfile=1:num_of_files
   
    fname=['tracks_ptv_' num2str(nfile) '.cam3.dat'];
    while (exist(fname,'file')==0)
        nfile=nfile+1;
        fname=['tracks_ptv_' num2str(nfile) '.cam3.dat'];
    end
    
    fid=fopen(fname,'rb');
    fread(fid,5,'char');
    Ntracks=fread(fid,1,'int32');
    fps=fread(fid,1,'float32');
    exp_time=fread(fid,1,'int32');
    treshold=fread(fid,1,'uchar');
    max_interp=fread(fid,1,'char');
    min_len=fread(fid,1,'int16');
    
    j=0;
    for nn=1:Ntracks;
        nnframes=fread(fid,1,'int32');
        if(nnframes>=minlen)
            j=j+1;
            nframes(j)=nnframes;
            %disp(sprintf('# of frames in track %i of %i : %i\n',nn,Ntracks,nframes(j)));
            for k=1:nframes(j)
                if (k==1) 
                    frame0(j)=fread(fid,1,'int32');
                    frame(j,1)=frame0(j);
                else
                    frame(j,k)=fread(fid,1,'int32');
                end
                x(j,k)=fread(fid,1,'float32');
                y(j,k)=fread(fid,1,'float32');
                z(j,k)=fread(fid,1,'float32');
                sigma_x(j,k)=fread(fid,1,'float32');
                sigma_y(j,k)=fread(fid,1,'float32');
                intensity(j,k)=fread(fid,1,'uchar');
                fake(j,k)=fread(fid,1,'char');
            end
            if nargin>4
                xf=conv(x(j,:),ker);
                yf=conv(y(j,:),ker);
                zf=conv(z(j,:),ker);
                
                %vx=conv(x(j,:),kerv);
                %vy=conv(y(j,:),kerv);
                %vz=conv(z(j,:),kerv);
                
                x(j,1:nframes(j)-2*filter_l)=xf(filter_l+1:nframes(j)-filter_l);
                x(j,nframes(j)-2*filter_l+1:nframes(j))=0;
                y(j,1:nframes(j)-2*filter_l)=yf(filter_l+1:nframes(j)-filter_l);
                y(j,nframes(j)-2*filter_l+1:nframes(j))=0;
                z(j,1:nframes(j)-2*filter_l)=zf(filter_l+1:nframes(j)-filter_l);
                z(j,nframes(j)-2*filter_l+1:nframes(j))=0;
                
                %vx(j,1:nframes(j)-2*filter_l)=vx(filter_l+1:nframes(j)-filter_l);
                %vx(j,nframes(j)-2*filter_l+1:nframes(j))=0;
                %vy(j,1:nframes(j)-2*filter_l)=vy(filter_l+1:nframes(j)-filter_l);
                %vy(j,nframes(j)-2*filter_l+1:nframes(j))=0;
                %vz(j,1:nframes(j)-2*filter_l)=vz(filter_l+1:nframes(j)-filter_l);
                %vz(j,nframes(j)-2*filter_l+1:nframes(j))=0;
                
                nframes(j)=nframes(j)-2*filter_l;
            end            
        else
            fseek(fid,nnframes*26,0);
        end
    end
    fclose(fid);
    
    if(j>0)
        Ntracks=j;
        lastframe=frame0+nframes-1;
        O0=ones(Ntracks,1)*lastframe-transpose(frame0)*ones(1,Ntracks)+1;
        [i,j]=find(O0.*O0'>0);
        O=zeros(Ntracks,Ntracks);
        for ii=1:size(i,1)
            O(i(ii),j(ii))=min([O0(i(ii),i(ii)),O0(j(ii),j(ii)),O0(i(ii),j(ii)),O0(j(ii),i(ii))]);
        end

        %O=O-diag(diag(O))
        O=triu(O,1);
        
        %2nd select tracks with sufficient long overlap
        [i_overlap,j_overlap]=find(O>minoverlap);
    
        if(~isempty(i_overlap))    
            for ii=1:size(i_overlap,1)
                len(ii_tot+ii)=O(i_overlap(ii),j_overlap(ii));
                nt(1:len(ii_tot+ii))=nt(1:len(ii_tot+ii))+1;
                if(frame0(i_overlap(ii))<frame0(j_overlap(ii)))
                    X1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=x(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
                    X2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=x(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
                    X1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    X2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Y1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=y(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
                    Y2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=y(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
                    Y1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Y2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Z1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=z(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
                    Z2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=z(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
                    Z1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Z2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                     if(nargin>4)
%                         VX1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vx(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
%                         VX2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vx(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
%                         VX1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VX2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VY1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vy(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
%                         VY2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vy(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
%                         VY1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VY2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VZ1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vz(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
%                         VZ2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vz(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
%                         VZ1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VZ2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                     end
                    %Y1(ii_tot+ii,:)=y(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+minoverlap);
                    %Y2(ii_tot+ii,:)=y(j_overlap(ii),1:minoverlap);                       %Z1(ii_tot+ii,:)=z(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+minoverlap);
                    %Z2(ii_tot+ii,:)=z(j_overlap(ii),1:minoverlap);
                    %pos(ii_tot+ii).X1(1:O(i_overlap(ii),j_overlap(ii)))=x(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
                    %pos(ii_tot+ii).X2(1:O(i_overlap(ii),j_overlap(ii)))=x(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));   
                else
                    X1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=x(i_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
                    X2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=x(j_overlap(ii),frame0(i_overlap(ii))-frame0(j_overlap(ii))+1:frame0(i_overlap(ii))-frame0(j_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
                    X1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    X2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Y1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=y(i_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
                    Y2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=y(j_overlap(ii),frame0(i_overlap(ii))-frame0(j_overlap(ii))+1:frame0(i_overlap(ii))-frame0(j_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
                    Y1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Y2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Z1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=z(i_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
                    Z2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=z(j_overlap(ii),frame0(i_overlap(ii))-frame0(j_overlap(ii))+1:frame0(i_overlap(ii))-frame0(j_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
                    Z1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
                    Z2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                     if(nargin>4)
%                         VX2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vx(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
%                         VX1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vx(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
%                         VX2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VX1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VY2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vy(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
%                         VY1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vy(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
%                         VY2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VY1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VZ2(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vz(i_overlap(ii),frame0(j_overlap(ii))-frame0(i_overlap(ii))+1:frame0(j_overlap(ii))-frame0(i_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));
%                         VZ1(ii_tot+ii,1:O(i_overlap(ii),j_overlap(ii)))=vz(j_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
%                         VZ2(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                         VZ1(ii_tot+ii,O(i_overlap(ii),j_overlap(ii))+1:maxoverlap)=0;
%                     end
                   %Y1(ii_tot+ii,:)=y(i_overlap(ii),1:minoverlap);
                   %Y2(ii_tot+ii,:)=y(j_overlap(ii),frame0(i_overlap(ii))-frame0(j_overlap(ii))+1:frame0(i_overlap(ii))-frame0(j_overlap(ii))+minoverlap);
                   %Z1(ii_tot+ii,:)=z(i_overlap(ii),1:minoverlap);
                   %Z2(ii_tot+ii,:)=z(j_overlap(ii),frame0(i_overlap(ii))-frame0(j_overlap(ii))+1:frame0(i_overlap(ii))-frame0(j_overlap(ii))+minoverlap);
                   %pos(ii_tot+ii).X1(1:O(i_overlap(ii),j_overlap(ii)))=x(i_overlap(ii),1:O(i_overlap(ii),j_overlap(ii)));
                   %pos(ii_tot+ii).X2(1:O(i_overlap(ii),j_overlap(ii)))=x(j_overlap(ii),frame0(i_overlap(ii))-frame0(j_overlap(ii))+1:frame0(i_overlap(ii))-frame0(j_overlap(ii))+O(i_overlap(ii),j_overlap(ii)));  
                end
                d2(ii_tot+ii,:)=(X1(ii_tot+ii,:)-X2(ii_tot+ii,:)).^2+(Y1(ii_tot+ii,:)-Y2(ii_tot+ii,:)).^2+(Z1(ii_tot+ii,:)-Z2(ii_tot+ii,:)).^2;
            end
            disp(sprintf('track file #%d done ; %d tracks found with %d overlapping points',nfile,size(i_overlap,1),minoverlap));
            ii_tot=ii_tot+ii;
        end
        %if(nargin>4)
            %eval(['save ' filename_save ' X1 X2 Y1 Y2 Z1 Z2 VX1 VY1 VZ1 VX2 VY2 VZ2 d2 len nt ;']);
        %else 
        if and(mod(nfile,5)==0,~isempty(i_overlap))
                eval(['save ' filename_save ' X1 X2 Y1 Y2 Z1 Z2 d2 len nt ;']);
        end
        clear x y z vx vy vz xf yf zf sigma_x sigma_y intensity fake frame0 lastframe nframes frame Ntracks nnframes O0 O i j i_overlap j_overlap;

        %end
    end
end