function R=xcorr_velocities(vel,bubble,nu0,theta)

c=340;
stat=stat_velocities(vel,nu0,theta);
fact=c/2/nu0/sin(theta/2)*32768;
[B,A]=butter(4,100*2/32768,'low');
dd=[];
for jj=1:numel(bubble.data)
    ddf=filtfilt(B,A,abs(bubble.data(jj).seg));
    bubble.data(jj).abs=abs(bubble.data(jj).seg)-ddf;
    dd=[dd ; bubble.data(jj).abs];
end;

mbub=mean(dd)
sbub=std(dd)

for ii=1:max(vel.length(vel.good))-1

    R(ii)=0;
    Nii=0;
    
    
    ind_j=find(vel.length(vel.good)>ii);
    
    N0=sum(vel.length(vel.good(ind_j)));
    
    %disp(sprintf('%i/%i : %i',ii,max(vel.length(vel.good))-1,numel(ind_j)));
    for jj=1:numel(ind_j)
        
        data=vel.data(vel.good(ind_j(jj))).freq;
        data=data*fact-stat.mean;
        
        Nii=Nii+length(data)-ii;
        datam=data(1:length(data)-ii);
        %datam=data(1+ii:length(data));
        %datap=abs(bubble.data(ind_j(jj)).seg(1+ii:length(data)))';
        datap=bubble.data(ind_j(jj)).abs(1+ii:length(data))-mbub;
        %datap=bubble.data(ind_j(jj)).abs(1:length(data)-ii)-mbub;
        %datap=data(1+ii:length(data));
                
        R(ii)=R(ii)+sum(datam.*datap');
        
    end
    R(ii)=R(ii)/Nii;
    %R(ii)=R(ii)*(N0-ii)/N0;
end
R=R/stat.std/sbub;