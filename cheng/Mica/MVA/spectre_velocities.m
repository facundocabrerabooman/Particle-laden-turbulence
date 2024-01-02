function [sp,f]=spectre_velocities(vel,varargin)

% [sp,f]=spectre_velocities(vel,varargin)
%
% computes power spectrum of the velocity signal in structure vel (field
% vel.data.freq
%
% Optional input parameters :
%
% varargin{1} = Sampling frequency
% varargin{2} = if set to 1 then computes only 1 fft for each "bubble"
%



if nargin>1
    Fs=varargin{1};
else
    Fs=1;
end

if nargin > 2
    vel.lmin=varargin{2};
else
    vel.lmin=1024;
end

if nargin>3
    fft_flag=varargin{3};
else
    fft_flag=0;
end

if nargin>4
    Ns=varargin{4};
else 
    Ns=vel.lmin;
end
if nargin>5
     stat=varargin{5};
%     c=340;
%     nu0=80000;
%     theta=165*pi/180;
%     stat=stat_velocities(vel,nu0,theta);
%     fact=c/2/nu0/sin(theta/2)*32768;
end

sp=zeros(Ns/2+1,1);

%ind=find((vel.length(vel.good)-150*3-1)>vel.lmin);
ll=zeros(numel(vel.good),1);
for j=1:numel(vel.good)
    ll(j)=length(vel.data(vel.good(j)).velf);
end

ind=find(ll>vel.lmin);
%ind=find(vel.length(vel.good)>vel.lmin);
%ind=find(vel.length(vel.good)>1024);
%ind=find(vel.length>vel.lmin);
%for j=1:length(vel.good)

h=spectrum.welch('Hamming',Ns,25,'UserDefined');
sp=zeros(Ns/2+1,1);

 for j=1:numel(ind)
       if fft_flag==1
         sig=vel.data(vel.good(ind(j))).velf(floor(vel.length(vel.good(ind(j)))/2)-floor(vel.lmin/2)+1:floor(vel.length(vel.good(ind(j)))/2)+floor(vel.lmin/2));
         
         if nargin>5
             m=stat.velf.mean;
         else
             m=mean(sig);
         end
             
         [p,f]=pwelch(sig-m,hanning(vel.lmin),[],vel.lmin,Fs);
         k(j)=1;
     else
         %sig=vel.data(vel.good(ind(j))).freq;
         sig=vel.data(vel.good(ind(j))).velf;
         
         if nargin>5
             m=stat.velf.mean;
         else
             m=mean(sig);
         end
        [p,f]=pwelch(sig-m,hanning(Ns),floor(Ns/4),Ns,Fs);
        %[p,f]=periodogram(sig-m,hanning(ll(ind(j))),Ns,Fs);
        %[p,f]=pmusic(sig-m,4,Ns,Fs);
        %p=psd(h,sig-m,'Fs',Fs,'NFFT',Ns);
        %[p,f]=spectrum(sig-m,Ns,floor(Ns/4),hanning(Ns),Fs);
       %[p,f]=spectrum(sig,vel.lmin,vel.lmin/4,hanning(vel.lmin),Fs);
       k(j)=fix((ll(ind(j))-floor(Ns/4))/(Ns-floor(Ns/4)));
     end
     
    %[p,f]=spectrum(vel.data(vel.good(ind(j))).freq,vel.lmin,vel.lmin/4,hanning(vel.lmin),Fs);
    %[p,f]=spectrum(vel.data(ind(j)).freq,vel.lmin,vel.lmin/4,blackman(vel.lmin),Fs);
    
    %[p,f]=pwelch(sig,hanning(vel.lmin),vel.lmin/4,vel.lmin,Fs);
    %[p,f]=pwelch(vel.data(vel.good(ind(j))).velf,hanning(vel.lmin),vel.lmin/4,vel.lmin,Fs);
    
    %sp=sp+p.data;
    sp=sp+p*k(j);
    loglog(f,sp);
end
 sp=sp/sum(k);
 %f=p.frequencies;