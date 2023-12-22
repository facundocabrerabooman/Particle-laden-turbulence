function [sp,f]=vel_spectrum(vel,Fsamp)
sp=zeros(129,1);
for j=1:length(vel)
    [p,f]=spectrum(vel(j).freq,256,64,hanning(256),Fsamp);
    sp=sp+p(:,1);
    %sp(j).f=f;
    %sp(j).sp=p(:,1);
end
sp=sp/length(vel);