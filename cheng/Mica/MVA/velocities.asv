function [vel]=velocities(fname,thresh_Hes,varargin)

%thresh_Hes

NbMic=7;
K=13;

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

if ~exist('vel')
    vel=struct();
    vel.good=[];
    vel.lmin=length_tresh;
end

for j=1:length(segments);
%for j=1:5
    %disp(sprintf('\n analyse de la bulle #%i',j));
    [freq,iHes]=MVA13mult1(conj(segments(j).seg'),1,NbMic,K);
    
    %vel(j+j_append).freq=freq(1,find(freq(1,:)>0))';
    %vel(j+j_append).Hes=iHes(find(freq(1,:)>0))';
    size(freq);
    freq=(freq(1,:));
    iHes=(iHes(1,:));
    ind=find(iHes<thresh_Hes);
    if(max(diff(ind))>1)
        indl=[find(diff(ind)>1)];
        if(iHes(ind(numel(ind)))<thresh_Hes)
            indl=[indl numel(ind)];
        end
        if (length(indl)==1)
            if(indl<length(ind)/2)
                m=length(ind)-indl;
                indl=[indl,length(ind)];
                im=1;
            else
                m=indl;
                indl=[0,indl];
                im=1;
            end
        else
            indl=[0 indl];
            l=diff(indl);
            [m,im]=max(l);
        end
    else
        m=length(ind);
        indl=[0,length(ind)];
        im=1;
    end
    %disp(sprintf('\n analyse de la bulle #%i ; m = %i',j,m));
  
    if(m > length_tresh)
        %disp(sprintf('\n analyse de la bulle #%i ; m = %i',j,m));
        vel.good=[vel.good j+j_append];
        ii=ind(indl(im)+1:indl(im+1));
        vel.data(j+j_append).seg=segments(j).seg(ii);
        vel.data(j+j_append).freq=freq(ii);
        vel.data(j+j_append).iHes=iHes(ii);
        vel.i0(j+j_append)=i0(j)+ii(1);
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

