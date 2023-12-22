function [tau,corra,corrv,arms,vrms,nc]=pump_magfield_corra_dt2(track,ndt,fps)

% [ca,tau,cv,arms,vrms,nc]=pump_magfield_corra_dt(res,ndt,ifplot);
% res is the structure containing data res(k).X res(k).Y res(k).Z
% should not contain nans
% ndt number of points over which dx and d2x are computed
% ifplot=0 for no plot
% vrms(i) rms of di with = 1,2,3 for x,y,z
% arms(i) rms of d2i 
% ca(i,:) correlation of d2i 
% cv(i,:) correlation of di 
%
% r.volk 2023/09/13

Lt=length(track);
% fps=1;%4000 if needed

xm=zeros(1,Lt);
ym=zeros(1,Lt);
zm=zeros(1,Lt);

vxm=zeros(1,Lt);
vym=zeros(1,Lt);
vzm=zeros(1,Lt);

vx2m=zeros(1,Lt);
vy2m=zeros(1,Lt);
vz2m=zeros(1,Lt);

axm=zeros(1,Lt);
aym=zeros(1,Lt);
azm=zeros(1,Lt);

ax2m=zeros(1,Lt);
ay2m=zeros(1,Lt);
az2m=zeros(1,Lt);

Lx=zeros(1,Lt);
Ly=zeros(1,Lt);
Lz=zeros(1,Lt);

for ktraj=1:Lt
    
    xx=track(ktraj).X;
    yy=track(ktraj).Y;
    zz=track(ktraj).Z;
    ll=length(zz);
    
    % ii=[1+ndt:ll-ndt];
                
    % be aware vx is not a velocity, should be devided by ndt
    vx=(xx(1+2*ndt(1):ll)-xx(1:end-2*ndt(1)))/2;
    vy=(yy(1+2*ndt(2):ll)-yy(1:end-2*ndt(2)))/2;
    vz=(zz(1+2*ndt(3):ll)-zz(1:end-2*ndt(3)))/2;
    
    % be aware ax is not an acceleration, should be devided by ndt^2
    ax=xx(1+2*ndt(1):ll)+xx(1:end-2*ndt(1))-2*xx(1+ndt(1):end-ndt(1));
    ay=yy(1+2*ndt(2):ll)+yy(1:end-2*ndt(2))-2*yy(1+ndt(2):end-ndt(2));
    az=zz(1+2*ndt(3):ll)+zz(1:end-2*ndt(3))-2*zz(1+ndt(3):end-ndt(3));
               
    Lx(ktraj)=length(vx);
    Ly(ktraj)=length(vy);
    Lz(ktraj)=length(vz);

    % mean for each track
    xm(ktraj)=mean(xx,"omitnan");
    ym(ktraj)=mean(yy,"omitnan");
    zm(ktraj)=mean(zz,"omitnan");
    
    vxm(ktraj)=mean(vx,"omitnan");
    vym(ktraj)=mean(vy,"omitnan");
    vzm(ktraj)=mean(vz,"omitnan");
    
    vx2m(ktraj)=mean(vx.^2,"omitnan");
    vy2m(ktraj)=mean(vy.^2,"omitnan");
    vz2m(ktraj)=mean(vz.^2,"omitnan");
    
    axm(ktraj)=mean(ax,"omitnan");
    aym(ktraj)=mean(ay,"omitnan");
    azm(ktraj)=mean(az,"omitnan");
    
    ax2m(ktraj)=mean(ax.^2,"omitnan");
    ay2m(ktraj)=mean(ay.^2,"omitnan");
    az2m(ktraj)=mean(az.^2,"omitnan");
    
    clear xx yy zz;
    
end

disp('computed moments');

Lmax=max([max(Lx),max(Ly),max(Lz)]);

% global mean
meanax=sum(Lx.*axm,"omitnan")/sum(Lx);
meanax2=sum(Lx.*ax2m,"omitnan")/sum(Lx);
meanay=sum(Ly.*aym,"omitnan")/sum(Ly);
meanay2=sum(Ly.*ay2m,"omitnan")/sum(Ly);
meanaz=sum(Lz.*azm,"omitnan")/sum(Lz);
meanaz2=sum(Lz.*az2m,"omitnan")/sum(Lz);

mx=sum(Lx.*xm,"omitnan")/sum(Lx);
my=sum(Ly.*ym,"omitnan")/sum(Ly);
mz=sum(Lz.*zm,"omitnan")/sum(Lz);

mvx=sum(Lx.*vxm,"omitnan")/sum(Lx);
mvx2=sum(Lx.*vx2m,"omitnan")/sum(Lx);
mvy=sum(Ly.*vym,"omitnan")/sum(Ly);
mvy2=sum(Ly.*vy2m,"omitnan")/sum(Ly);
mvz=sum(Lz.*vzm,"omitnan")/sum(Lz);
mvz2=sum(Lz.*vz2m,"omitnan")/sum(Lz);

axrms=sqrt(meanax2-meanax^2);
vxrms=sqrt(mvx2-mvx^2);
ayrms=sqrt(meanay2-meanay^2);
vyrms=sqrt(mvy2-mvy^2);
azrms=sqrt(meanaz2-meanaz^2);
vzrms=sqrt(mvz2-mvz^2);

arms.x = axrms;
arms.y = ayrms;
arms.z = azrms;
vrms.x = vxrms;
vrms.y = vyrms;
vrms.z = vzrms;

cax=zeros(1,2*Lmax-1);
cay=zeros(1,2*Lmax-1);
caz=zeros(1,2*Lmax-1);

cvx=zeros(1,2*Lmax-1);
cvy=zeros(1,2*Lmax-1);
cvz=zeros(1,2*Lmax-1);

nc.x=zeros(size(cax));
nc.y=zeros(size(cay));
nc.z=zeros(size(caz));

m=0;

for ktraj=1:length(track)
    
    xx=track(ktraj).X;
    yy=track(ktraj).Y;
    zz=track(ktraj).Z;
    ll=length(zz);
    
    % ii=[1+ndt:ll-ndt];    
        
    % be aware vx is not a velocity, should be devided by ndt
    vx=(xx(1+2*ndt(1):ll)-xx(1:end-2*ndt(1)))/2;
    vy=(yy(1+2*ndt(2):ll)-yy(1:end-2*ndt(2)))/2;
    vz=(zz(1+2*ndt(3):ll)-zz(1:end-2*ndt(3)))/2;
    
    % be aware ax is not an acceleration, should be devided by ndt^2
    ax=xx(1+2*ndt(1):ll)+xx(1:end-2*ndt(1))-2*xx(1+ndt(1):end-ndt(1));
    ay=yy(1+2*ndt(2):ll)+yy(1:end-2*ndt(2))-2*yy(1+ndt(2):end-ndt(2));
    az=zz(1+2*ndt(3):ll)+zz(1:end-2*ndt(3))-2*zz(1+ndt(3):end-ndt(3));

    vx2=vx'-mvx;
    vy2=vy'-mvy;
    vz2=vz'-mvz;
    clear vx vy vz
    
    ax2=ax'-meanax;
    ay2=ay'-meanay;
    az2=az'-meanaz;
    
    jjx=length(ax2);
    jjy=length(ay2);
    jjz=length(az2);
    
    % computes distance to the mean position
%     dist=mean(sqrt((xx-mx).^2+(yy-my).^2+(zz-mz).^2));

    c=xcorr(ax2,'none');
    cax(Lmax-jjx+1:Lmax+jjx-1)=cax(Lmax-jjx+1:Lmax+jjx-1)+c;
    c=xcorr(ay2,'none');
    cay(Lmax-jjy+1:Lmax+jjy-1)=cay(Lmax-jjy+1:Lmax+jjy-1)+c;
    c=xcorr(az2,'none');
    caz(Lmax-jjz+1:Lmax+jjz-1)=caz(Lmax-jjz+1:Lmax+jjz-1)+c;
    
    c=xcorr(vx2,'none');
    cvx(Lmax-jjx+1:Lmax+jjx-1)=cvx(Lmax-jjx+1:Lmax+jjx-1)+c;
    c=xcorr(vy2,'none');
    cvy(Lmax-jjy+1:Lmax+jjy-1)=cvy(Lmax-jjy+1:Lmax+jjy-1)+c;
    c=xcorr(vz2,'none');
    cvz(Lmax-jjz+1:Lmax+jjz-1)=cvz(Lmax-jjz+1:Lmax+jjz-1)+c;
    
    clear c;
    
    nc.x(Lmax-jjx+1:Lmax+jjx-1)=nc.x(Lmax-jjx+1:Lmax+jjx-1)+[1:jjx,jjx-1:-1:1];
    nc.y(Lmax-jjy+1:Lmax+jjy-1)=nc.y(Lmax-jjy+1:Lmax+jjy-1)+[1:jjy,jjy-1:-1:1];
    nc.z(Lmax-jjz+1:Lmax+jjz-1)=nc.z(Lmax-jjz+1:Lmax+jjz-1)+[1:jjz,jjz-1:-1:1];
    m=m+1;
    
    clear xx yy zz ax ay az vx vy vz;
    
    if round(ktraj/5000)*5000==ktraj
       disp(ktraj);
    end
       
end        

t=(-Lmax+1:Lmax-1)/fps;
Jx=find(nc.x~=0);
Jy=find(nc.y~=0);
Jz=find(nc.z~=0);
corrax=cax(Jx)./nc.x(Jx);
corrvx=cvx(Jx)./nc.x(Jx);
corray=cay(Jy)./nc.y(Jy);
corrvy=cvy(Jy)./nc.y(Jy);
corraz=caz(Jz)./nc.z(Jz);
corrvz=cvz(Jz)./nc.z(Jz);

taux=t(Jx);
tauy=t(Jy);
tauz=t(Jz);

taux1=taux((length(taux)-1)/2+1:end);
tauy1=tauy((length(tauy)-1)/2+1:end);
tauz1=tauz((length(tauz)-1)/2+1:end);
corrax1=corrax((length(corrax)-1)/2+1:end);
corrvx1=corrvx((length(corrvx)-1)/2+1:end);
corray1=corray((length(corray)-1)/2+1:end);
corrvy1=corrvy((length(corrvy)-1)/2+1:end);
corraz1=corraz((length(corraz)-1)/2+1:end);
corrvz1=corrvz((length(corrvz)-1)/2+1:end);


clear taux corrax corrvx
clear tauy corray corrvy
clear tauz corraz corrvz

corra.x = corrax1;
corra.y = corray1;
corra.z = corraz1;

corrv.x = corrvx1;
corrv.y = corrvy1;
corrv.z = corrvz1;

tau.x = taux1;
tau.y = tauy1;
tau.z = tauz1;


%disp([int2str(m) ' trajectories in a ball of size ' num2str(ballsize)])

