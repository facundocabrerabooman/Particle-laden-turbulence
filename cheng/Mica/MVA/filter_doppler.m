function [fdata,FDoppler,WDoppler,SpDoppler]=filter_doppler(data,Fsamp)

%[fdata,FDoppler,WDoppler]=filter_doppler(data,F_Ech);

% st=fft(data);
% st(round(length(data)/2)+1:length(data))=0;
% data=ifft(st);
% %[sp,f]=spectrum(real(data),4096,1024,hanning(4096),Fsamp);
[sp,f]=pwelch((data-mean(data)),hanning(4096),1024,4096,Fsamp,'centered');
%[sp_complex,f]=pwelch(real(data-mean(data)),hanning(4096),1024,4096,Fsamp);
% [B,A]=butter(4,0.25,'low');
%sp=filtfilt(B,A,sp);
sp=medfilt1(sp,10);
 %figure;semilogy(f,sp(:,1));
%return;
II=find(f<50);
f(II)=[];
sp(II)=[];
%sp_complex(II)=[];

pk=1;
n=1;
while(pk<50)
    [m,pk]=max(sp(n:length(sp(:,1)),1));
    n=n+1;
end
pk=pk+n-1;

[m,im]=min(sp(1:pk,1));
width=2*(pk-im)*Fsamp/4096;
FDoppler=f(pk);
SpDoppler=sp(pk);
Fmin=FDoppler-width;
Fmax=FDoppler+3*width;
WDoppler=Fmax-Fmin;

 [B,A]=butter(4,[max(1.1*min(f),Fmin) min(Fmax,.9*max(f))]*2/Fsamp,'bandpass');
 %[B,A]=butter(4,[7000 9000]*2/Fsamp,'bandpass');
%[B,A]=cheby1(2,0.5,[FDoppler-WDoppler/2 FDoppler+WDoppler/2]*2/Fsamp);
fdata=filtfilt(B,A,data);
%st=fft(fdata);
%st(round(length(fdata)/2)+1:length(fdata))=0;
%fdata=ifft(st);