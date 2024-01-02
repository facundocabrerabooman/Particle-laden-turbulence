function [vel]=velocities_bubble(bubble,varargin)


length_tresh=1024;
j_append=0;
if nargin>1
    vel=varargin{1};
    j_append=length(vel.data);
end

%[data,F_Span,Fsamp,F_Cent,Range,Scale_Factor,Ovld_Flag]=lire_e1430_cmplx(fname,1,1024*1024,0);

%[fdata,FDoppler,WDoppler]=filter_doppler(data,Fsamp);

%fdata=decimate(fdata,2);
%[segments,i0,id]=get_seg(fdata);


if ~exist('vel')
    vel=struct();
    vel.good=[];
    vel.lmin=length_tresh;
end

for j=1:numel(bubble.length);
%for j=1:5
    %disp(sprintf('\n analyse de la bulle #%i',j));
    [freq,iHes]=MVA13mult1(flipud(fliplr(bubble.data(j).fseg')),1);
    %vel(j+j_append).freq=freq(1,find(freq(1,:)>0))';
    %vel(j+j_append).Hes=iHes(find(freq(1,:)>0))';
    size(freq);
    freq=fliplr(freq(1,:));
    iHes=fliplr(iHes(1,:));
    ind=find(iHes<1.e-5);
    indl=find(diff(ind)>1);
    l=diff(indl);
    [m,im]=max(l);
    if(m > length_tresh)
        vel.good=[vel.good j+j_append];
        ii=ind(indl(im)+1:indl(im+1));  
        vel.data(j+j_append).freq=freq(ii);
        vel.data(j+j_append).iHes=iHes(ii);
        vel.i0(j+j_append)=1+ind(1);
        vel.length(j+j_append)=m;
        %vel.file(j+j_append)=cellstr(fname);
        vel.status(j+j_append)=1;
    else
        vel.data(j+j_append).freq=freq;
        vel.data(j+j_append).iHes=iHes;
        vel.i0(j+j_append)=1;
        vel.length(j+j_append)=bubble.length(j);
        %vel.file(j+j_append)=cellstr(fname);
        vel.status(j+j_append)=0;
    end
end

disp(sprintf('**** %i bonnes bulles \n',length(vel.good)));

%[sp,f]=vel_spectrum(vel,Fsamp);

