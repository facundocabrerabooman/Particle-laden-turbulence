function [vel]=velocities(fname,thresh_Hes,varargin)
thresh_Hes
length_tresh=512;
j_append=0;
if nargin>2
    vel=varargin{1};
    j_append=length(vel.data);
end

[data,F_Span,Fsamp,F_Cent,Range,Scale_Factor,Ovld_Flag]=lire_e1430_cmplx(fname,1,1024*1024,0);

[fdata,FDoppler,WDoppler]=filter_doppler(data,Fsamp);

%fdata=decimate(fdata,2);
[segments,i0,id]=get_seg(fdata);

if isempty(segments)
    return;
end

if ~exist('vel')
    vel=struct();
    vel.good=[];
    vel.lmin=length_tresh;
end

% [p,f]=pwelch(fdata,hanning(length(fdata)));
% ff=f-pi;
% pp=fftshift(p);
% tt=max(pp(find(ff<0)));
% sigeps2=tt/2;
sigeps2=1e-10;
    
for j=1:length(segments);
%for j=1:5

    %estimation de sigeps2;
    %[p,f]=pwelch(segments(j))%.seg,hanning(1024),512,1024,1);
    %figure;plot(f-0.5,fftshift(p));
%     ff=f-0.5;
%     pp=fftshift(p);
%     tt=max(pp(find(ff<0)));
%     sigeps2=tt/2;
    %disp(sprintf('\n analyse de la bulle #%i',j));
    [freq,iHes]=MVA13mult1(conj(segments(j).seg'),sigeps2,1);
    %vel(j+j_append).freq=freq(1,find(freq(1,:)>0))';
    %vel(j+j_append).Hes=iHes(find(freq(1,:)>0))';
    size(freq);
    freq=(freq(1,:));
    iHes=(iHes(1,:));
    ind=find(iHes<thresh_Hes);
    if(max(diff(ind))>1)
        indl=find(diff(ind)>1);
        l=diff(indl);
        [m,im]=max(l);
    else
        m=length(ind);
        indl=[ind(1),ind(length(ind))-1];
        im=1;
    end
    
    if(m > length_tresh)
        vel.good=[vel.good j+j_append];
        ii=ind(indl(im)+1:indl(im+1));
        vel.data(j+j_append).seg=segments(j).seg(ii);
        vel.data(j+j_append).freq=freq(ii);
        vel.data(j+j_append).iHes=iHes(ii);
        vel.i0(j+j_append)=i0(j)+ind(1);
        %vel.length(j+j_append)=m;
        vel.length(j+j_append)=length(freq(ii));
        vel.file(j+j_append)=cellstr(fname);
        vel.status(j+j_append)=1;
    else
        vel.data(j+j_append).seg=segments(j).seg;
        vel.data(j+j_append).freq=freq;
        vel.data(j+j_append).iHes=iHes;
        vel.i0(j+j_append)=i0(j);
        %vel.length(j+j_append)=id(j);
        vel.length(j+j_append)=length(freq);
        vel.file(j+j_append)=cellstr(fname);
        vel.status(j+j_append)=0;
    end
end

disp(sprintf('**** %i bonnes bulles \n',length(vel.good)));

%[sp,f]=vel_spectrum(vel,Fsamp);

