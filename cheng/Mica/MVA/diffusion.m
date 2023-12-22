function [D,N]=diffusion(vel,stat)

%c=340;
%fact=c/2/nu0/sin(theta/2)*32768;
%fm=stat.velf.mean/fact;

for j=1:numel(vel.good)
    len(j)=length(vel.data(vel.good(j)).velf);
end

t=linspace(0,(max(len)-1)/32768,max(len));
N=zeros(1,max(len));
D=zeros(numel(vel.good,max(len)));

for j=1:numel(vel.good)
    for k=1:len(j)
        D(j,k)=integ(detrend(vel.data(vel.good(j)).velf(1:k)-stat.velf.mean),t(1:k));
    end
    N(1:len(j))=N(1:len(j))+1;
end