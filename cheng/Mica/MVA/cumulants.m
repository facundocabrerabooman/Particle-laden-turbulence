function [M,C]=cumulants(vel_in,p,dt)

%[M]=cumulants(vel,p,dt)
%
% calcule les moments d'ordre p des incrŽments 
%

M=zeros(numel(p),numel(dt));

for kk=1:numel(dt)
    lc=0;
    for jj=1:numel(vel_in.good)
        if(length(vel_in.data(vel_in.good(jj)).velf)-dt(kk)>1)
            ind=1:length(vel_in.data(vel_in.good(jj)).velf)-dt(kk);
            x=vel_in.data(vel_in.good(jj)).velf(ind+dt(kk))-vel_in.data(vel_in.good(jj)).velf(ind);
            
            for ii=1:numel(p)
                M(ii,kk)=M(ii,kk)+sum(abs(x).^p(ii));
            end
                
            %data(ig).l=numel(data(ig).x);
            lc=lc+numel(x);
            %data(ig).lc=lc;
            %inc=[inc vel_in.data(vel_in.good(jj)).velf(ind+dt(kk))-vel_in.data(vel_in.good(jj)).velf(ind)];
        end
    end
    M(:,kk)=M(:,kk)/lc;
    
    P=polyfit(p,log(M(:,kk))'./p,2);
    C.c1(kk)=P(3);
    C.c2(kk)=P(2);
    C.c3(kk)=P(1);
end