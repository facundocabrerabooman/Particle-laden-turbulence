function sp=sp_sawford(f,TL,Re,stat);


% sp=sp_sawford(f,TL,Re,stat);

%c=340;
%nu0=80000;
%theta=165*pi/180;
%stat=stat_velocities(vel,nu0,theta);
%fact=c/2/nu0/sin(theta/2)*32768;


T2=TL/sqrt(Re);
sp=2*(stat.velf.std)^2*(TL^2./(1+(2*pi*f*TL).^2)-T2^2./(1+(2*pi*f*T2).^2))/(TL-T2);