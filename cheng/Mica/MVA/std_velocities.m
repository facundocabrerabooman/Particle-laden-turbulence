function s=std_velocities(vel,stat)

%c=340;
%stat=stat_velocities(vel,nu0,theta);
%fact=c/2/nu0/sin(theta/2)*32768;
length_velf=zeros(1,numel(vel.data));
for kk=1:numel(vel.good)
    length_velf(vel.good(kk))=length(vel.data(kk).velf);
end

s=0.;
N=0;

for jj=1:numel(vel.good)
        data=[];
        datam=[];
        datap=[];
        data=vel.data(vel.good(jj)).velf;
        data=data-stat.velf.mean;
        %datas=[datas data];
        
        N=N+length(data);
                
        s=s+sum(data.*data);
   
end
    
    s=s/(N-1);
