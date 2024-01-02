function vel_out=acc_pente(vel,F_Ech,N)

 % taille de la fenetre de fit lineaire ; prendre N impair
%acc=zeros(30,1);
%acc=zeros(vel.length(vel.good(1))-(N-1),1);

vel_out=vel;
for jj=1:numel(vel.good);
    acc=[];
    disp(sprintf('%i',jj));
    for k=(N-1)/2+1:(vel.length(vel.good(jj)))-(N-1)/2;
        %           j=N/2+1:(vel.length(vel.good(1))-(N/2+1));
        ind=[k-(N-1)/2:k+(N-1)/2];
        s=vel.data(vel.good(jj)).freq(ind);
        P=polyfit(ind/F_Ech,s,1);
        acc=[acc P(1)];
        %acc(k-(N-1)/2)=P(1);
    end
    vel_out.data(vel.good(jj)).acc=acc;
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

% function vel_out=accelerations(vel_in,w)
% 
% % calculates acceleration from velocities in structure vel_in
% % the gaussian kernel is calculated from 
% % kernel = velfiltcoef(filterwidth, filterlen)
% 
% kerv=posfiltcoef(w,3*w);
% kera=velfiltcoef(w,3*w);
% 
% 
% %l=length(kera);
% vel_out=vel_in;
% for jj=1:numel(vel_in.data)
%    disp(sprintf('%i',jj));
%    
%    acctmp=conv(kera,vel_in.data(jj).freq);
%    %acc=acctmp(l:length(acctmp)-l);
%    
%    veltmp=conv(kerv,vel_in.data(jj).freq);
%    %velf=veltmp(l:length(veltmp)-l);
%    
%    vel_out.data(jj).acc=acctmp(numel(kera):vel_in.length(jj));
%    vel_out.data(jj).velf=veltmp(numel(kerv):vel_in.length(jj));
%    
%    %for kk=1:vel_in.length(jj)-numel(ker)
%    %     acc(kk)=sum(ker'.*vel_in.data(jj).freq(kk:kk+numel(ker)-1));
%    %     vel.data(jj).acc=acc;
%    %end
% end