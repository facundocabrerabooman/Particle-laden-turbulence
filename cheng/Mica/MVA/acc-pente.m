function vel_out=acc-pente(vel,Fc,F_Ech)

N=10:10:200;

for i=N/2+1:vel.length(vel.good(j)-(N/2+1))
    ind=[j-N/2,j+N/2];
    s=vel.data.vel.good(j).freq(ind);
    p=polyfit(ind/F_Ech,s,1);
    a(j-N/2)=p(1);
    vel.data(vel.good(j)).acc=a
end
    
    


% function [sigma_a, sigma_v, w]=findFilterWidth(vel);
% 
% w=5:5:250;
% l=3*w;
% 
% 
% for j=1:numel(w)
%     disp(sprintf('w = %i',w(j)));
%     kerp = posfiltcoef(w(j),l(j));
%     kerv = velfiltcoef(w(j),l(j));
%     velf=[];
%     acc=[];
%     for jj=1:numel(vel.good)
%         velftmp=conv(kerp,vel.data(vel.good(jj)).freq);
%         acctmp=conv(kerv,vel.data(vel.good(jj)).freq);
%         
%         velf=[velf velftmp(l(j):length(velftmp)-l(j))];
%         acc=[acc acctmp(l(j):length(acctmp)-l(j))];
%         clear velftmp acctmp;
%     end
%     
%     sigma_v(j)=std(velf);
%     sigma_a(j)=std(acc);
% end