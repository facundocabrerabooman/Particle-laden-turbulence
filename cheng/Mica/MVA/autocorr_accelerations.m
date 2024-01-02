function [R,S,N]=autocorr_accelerations(vel,stat)

%c=340;
%stat=stat_velocities(vel,nu0,theta);
%fact=c/2/nu0/sin(theta/2)*32768;

length_acc=zeros(1,numel(vel.data));
for kk=1:numel(vel.good)
    length_acc(vel.good(kk))=length(vel.data(vel.good(kk)).acc);
end

for ii=0:max(length_acc(vel.good))-1    
    S(ii+1)=0.;
    R(ii+1)=0.;
    N(ii+1)=0;
    %Sp(ii)=0.;
    %Sm(ii)=0.;
    %datas=[];
    
    ind_j=find(length_acc(vel.good)>ii);
    
    N0=sum(length_acc(vel.good(ind_j)));
    
    %disp(sprintf('%i/%i : %i',ii,max(vel.length(vel.good))-1,numel(ind_j)));
    for jj=1:numel(ind_j)
        data=[];
        datam=[];
        datap=[];
        data=vel.data(vel.good(ind_j(jj))).acc;
        data=data-stat.acc.mean;
        %datas=[datas data];
        
        N(ii+1)=N(ii+1)+length(data)-ii;
        datam=data(1:length(data)-ii);
        datap=data(1+ii:length(data));
                
        R(ii+1)=R(ii+1)+sum(datam.*datap);
        S(ii+1)=S(ii+1)+sum((datap-datam).^2);
        %Sp(ii)=Sp(ii)+sum(datap.^2);
        %Sm(ii)=Sm(ii)+sum(datam.^2);
    end
    
    R(ii+1)=R(ii+1)/N(ii+1);
    S(ii+1)=S(ii+1)/N(ii+1);
    %Sm(ii)=Sm(ii)/N(ii);
    %Sp(ii)=Sp(ii)/N(ii);
    %sigma(ii)=std(datas);
    %moy(ii)=mean(datas);
    
    %R(ii)=R(ii)*(N0-ii)/N0;
end
R=R/stat.acc.std^2;