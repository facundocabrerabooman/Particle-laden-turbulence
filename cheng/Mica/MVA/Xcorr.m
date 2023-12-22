function [R,S]=Xcorr(sig1,sig2)

%[R,S]=Xcorr(sig)
% R est calculee par R(kk+1)=mean(sig1(ind+kk).*sig2(ind));
% equivalent a xcorr avec l'option 'unbiased'

data1=[];
data2=[];
data1=detrend(sig1,'constant');
data2=detrend(sig2,'constant');
sig=std(data1).^2;
for kk=0:numel(data1)-1
	ind=1:numel(data1)-kk;
	R(kk+1)=mean(data2(ind+kk).*data1(ind));
	S(kk+1)=mean((data2(ind+kk)-data2(ind)).*(data1(ind+kk)-data1(ind)));
end
R=R/sig;