function [sp,f]=spectre_acoustAmp(bubble,varargin)

vel.lmin=1024;

if nargin>1
    Fs=varargin{1};
else
    Fs=1;
end


%lmin=min(bubble.length)
lmin=1024;
sp=zeros(floor(lmin/2)+1,1);
%ind=find(vel.length(vel.good)>vel.lmin);
%ind=find(vel.length>vel.lmin);
%for j=1:length(vel.good)
 for j=1:numel(bubble.length)
    [p,f]=spectrum(abs(bubble.data(j).seg),lmin,lmin/4,blackman(lmin),Fs);
    %[p,f]=spectrum(vel.data(ind(j)).freq,vel.lmin,vel.lmin/4,blackman(vel.lmin),Fs);
    sp=sp+p(:,1);
end
sp=sp/numel(bubble.length);