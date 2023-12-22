function res=bpass2(Iimg,lnoise,lobject,Outside_level)

% res=bpass2(Iimg,lnoise,lobject)

if ~exist('Outside_level','var')
    Outside_level=mean(mean(Iimg));
end
[npixy, npixx] = size(Iimg);
DD=2*lobject+1;

res=zeros(npixy+DD,npixx+DD)+Outside_level;
res((DD-1)/2+1:(DD-1)/2+npixy,(DD-1)/2+1:(DD-1)/2+npixx)=Iimg;

res=bpass(res,lnoise,lobject);

res=res((DD-1)/2+1:(DD-1)/2+npixy,(DD-1)/2+1:(DD-1)/2+npixx);
res=res*254/max(max(res));