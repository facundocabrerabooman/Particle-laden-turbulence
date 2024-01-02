function [vel]=velocities2(fname,varargin)

%thresh_Hes

NbMic=7;
K=13;

length_tresh=512;
j_append=0;
if nargin>1
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
    [freq,iHes]=MVA13mult1(conj(segments(j).seg'),1,NbMic,K);
    vel.data(j+j_append).seg=segments(j).seg;
    vel.data(j+j_append).freq=freq;
    vel.data(j+j_append).iHes=iHes;
    vel.i0(j+j_append)=i0(j);
    vel.length(j+j_append)=length(freq);
    vel.file(j+j_append)=cellstr(fname);
    vel.status(j+j_append)=0;
end
disp(sprintf('**** %i bulles \n',length(vel.data)));

%[sp,f]=vel_spectrum(vel,Fsamp);

