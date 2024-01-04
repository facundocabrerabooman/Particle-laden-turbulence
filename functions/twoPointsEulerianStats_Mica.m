function eulerStats= twoPointsEulerianStats_Mica(vtracks,dRminmax,Nbins,varargin)

%  Converts tracks from Clement matrices to vtrack structure
%  Colonne 1 : position x (mm)
%  2 : position y (mm)
%  3: numéro de l'image
%  4: numéro de la trajectoire
%  5 : rayon (mm)
%  6 : vx
%  7: vy
%  8: v=sqrt(vx^2+vy^2 )
%  9: theta (rad) sachant que le 0 correspond à une nage vers le bas de l'image
%  px   = 0.29 pour 1 5 10 20 30 50
%       = 0.23 pour 29 43

if nargin>4
    Nmax = varargin{2};
end

if nargin > 3
    px = varargin{1};
else
    px = 1;
end

dRmin = dRminmax(1);
dRmax = dRminmax(2);

%%
 w=10;
 L=3*w;
% kerp = posfiltcoef(w,L);
% kerv = velfiltcoef(w,L);
% kera = accfiltcoef(w,L);
% 
% if ~exist('Nmax')
%     vtracks = addStructFun(vtracks,'X','Xf',@(X)(conv(X,kerp,'valid')*px));
%     vtracks = addStructFun(vtracks,'Y','Yf',@(X)(conv(X,kerp,'valid')*px));
%     vtracks = addStructFun(vtracks,'X','Vx',@(X)(conv(X,kerv,'valid')*px*fps));
%     vtracks = addStructFun(vtracks,'Y','Vy',@(X)(conv(X,kerv,'valid')*px*fps));
%     vtracks = addStructFun(vtracks,'X','Ax',@(X)(conv(X,kera,'valid')*px*fps^2));
%     vtracks = addStructFun(vtracks,'Y','Ay',@(X)(conv(X,kera,'valid')*px*fps^2));
%     vtracks = addStructFun(vtracks,'T','Tf',@(X)(X(floor(L/2)+1:end-floor(L/2))));
% else
%     vtracks = addStructFun(vtracks,'X','Xf',@(X)(conv(X(1:Nmax),kerp,'valid')*px));
%     vtracks = addStructFun(vtracks,'Y','Yf',@(X)(conv(X(1:Nmax),kerp,'valid')*px));
%     vtracks = addStructFun(vtracks,'X','Vx',@(X)(conv(X(1:Nmax),kerv,'valid')*px*fps));
%     vtracks = addStructFun(vtracks,'Y','Vy',@(X)(conv(X(1:Nmax),kerv,'valid')*px*fps));
%     vtracks = addStructFun(vtracks,'X','Ax',@(X)(conv(X(1:Nmax),kera,'valid')*px*fps^2));
%     vtracks = addStructFun(vtracks,'Y','Ay',@(X)(conv(X(1:Nmax),kera,'valid')*px*fps^2));
%     vtracks = addStructFun(vtracks,'T','Tf',@(X)(X(floor(L/2)+1:Nmax-floor(L/2))));  
% end

eulerStats.filtW = w;
eulerStats.filtL = L;

%% Calculate acceleration

%Thresh =7;
%px = 0.23;%29 - 43
% px = 0.29;%50 - 30
%%
part = track2part(vtracks,{'Tf','Xf','Yf','Zf','Vx','Vy','Vz','Ax','Ay','Az'},1);

part = addStructFun(part,'Xf','X',@(X)(X));
part = addStructFun(part,'Yf','Y',@(X)(X));
part = addStructFun(part,'Zf','Z',@(X)(X));

Vx=reshape(vertcat(part.Vx),1,[]);
Vy=reshape(vertcat(part.Vy),1,[]);
Vz=reshape(vertcat(part.Vz),1,[]);
sigmaVx=std(Vx);
sigmaVy=std(Vy);
sigmaVz=std(Vz);
sigmaV=sqrt(sigmaVx^2+sigmaVy^2+sigmaVz.^2);

%%
Vmoy = arrayfun(@(X)(mean(sqrt(X.Vx.^2+X.Vy.^2+X.Vz.^2))),part);
VmoyX = arrayfun(@(X)(mean(X.Vx)),part);
VmoyY = arrayfun(@(X)(mean(X.Vy)),part);
VmoyZ = arrayfun(@(X)(mean(X.Vz)),part);
Vstd = arrayfun(@(X)(std(sqrt(X.Vx.^2+X.Vy.^2+X.Vz.^2))),part);
VstdX= arrayfun(@(X)(std(X.Vx)),part);
VstdY= arrayfun(@(X)(std(X.Vy)),part);
VstdZ= arrayfun(@(X)(std(X.Vz)),part);
Amoy = arrayfun(@(X)(mean(sqrt(X.Ax.^2+X.Ay.^2+X.Az.^2))),part);
AmoyX = arrayfun(@(X)(mean(X.Ax)),part);
AmoyY = arrayfun(@(X)(mean(X.Ay)),part);
AmoyZ = arrayfun(@(X)(mean(X.Az)),part);
Astd = arrayfun(@(X)(std(sqrt(X.Ax.^2+X.Ay.^2+X.Az.^2))),part);
AstdX= arrayfun(@(X)(std(X.Ax)),part);
AstdY= arrayfun(@(X)(std(X.Ay)),part);
AstdZ= arrayfun(@(X)(std(X.Az)),part);

eulerStats.Vmoy=Vmoy;
eulerStats.VmoyX=VmoyX;
eulerStats.VmoyY=VmoyY;
eulerStats.VmoyZ=VmoyZ;

eulerStats.Vstd=Vstd;
eulerStats.VstdX=VstdX;
eulerStats.VstdY=VstdY;
eulerStats.VstdZ=VstdZ;

eulerStats.Astd=Astd;
eulerStats.AstdX=AstdX;
eulerStats.AstdY=AstdY;
eulerStats.AstdZ=AstdZ;

%%
pairAll = pairStat3D(part);
pair = pairAll;
%% consider only portions of tracks where the average velocity Vmoy is
%% larger than threshold Thresh
%partSorted = sortPart(vtracks,Vmoy,.75);
%pairSorted = pairStat(partSorted);
%pair = pairSorted;
%%

%I=1:6286;

dX=[pair.dX];
dY=[pair.dY];
dZ=[pair.dZ];
dR2=[pair.dR2];
dVx=[pair.dVx];
dVy=[pair.dVy];
dVz=[pair.dVz];
dV2=[pair.dV2];
dAx=[pair.dAx];
dAy=[pair.dAy];
dAz=[pair.dAz];
dA2=[pair.dA2];

%% Calculate longitudinal increments
% [theta,r]=cart2pol(dX,dY);
% dA_trans=-sin(theta).*dAx+cos(theta).*dAy;
% dV_trans=-sin(theta).*dVx+cos(theta).*dVy;


dA_long=(dAx.*dX+dAy.*dY+dAz.*dZ)./(sqrt(dR2));
dV_long=(dVx.*dX+dVy.*dY+dVz.*dZ)./(sqrt(dR2));

deps = dAx.*dVx + dAy.*dVy + dAz.*dVz;
deps_long=dA_long.*dV_long;

%% Calculate structure functions
%dRmin = 5; %mm
%dRmax = 110;%mm
%Nbins = 40;

%[S2long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',dV_long.^2,'Fun',@mean);
%[S3x dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',dVx.^3,'Fun',@mean);
%[S3y dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',dVy.^3,'Fun',@mean);

% 
% [S2x dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dVx).^2,'Fun',@mean);
% [S2y dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dVy).^2,'Fun',@mean);
% 
% [S1long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',dV_long,'Fun',@mean);
% [S1longAbs dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dV_long),'Fun',@mean);
% [S2long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',dV_long.^2,'Fun',@mean);
% [S3long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',(dV_long).^3,'Fun',@mean);
% [S3longAbs dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dV_long).^3,'Fun',@mean);
% [S4long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',dV_long.^4,'Fun',@mean);
% [S5long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',(dV_long.^5),'Fun',@mean);
% [S5longAbs dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dV_long.^5),'Fun',@mean);
% 
% [Sau dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',deps,'Fun',@mean);
% [Sau_long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',deps_long,'Fun',@mean);
% 
% r=dR_bins{1};

dRmin = dRmin^2;
dRmax = dRmax^2;

[Nsamples dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',dVx,'Fun',@numel);

[S2x dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dVx).^2,'Fun',@mean);
[S2y dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dVy).^2,'Fun',@mean);
[S2z dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dVz).^2,'Fun',@mean);

[S1long dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',dV_long,'Fun',@mean);
[S1longAbs dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dV_long),'Fun',@mean);
[S2long dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',dV_long.^2,'Fun',@mean);
[S3long dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',(dV_long).^3,'Fun',@mean);
[S3longAbs dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dV_long).^3,'Fun',@mean);
[S4long dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',dV_long.^4,'Fun',@mean);
[S5long dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',(dV_long.^5),'Fun',@mean);
[S5longAbs dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',abs(dV_long.^5),'Fun',@mean);

[Sau dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',deps,'Fun',@mean);
[Sau_long dR_edges dR_bins]=histcn((dR2'),logspace(log10(dRmin),log10(dRmax),Nbins),'AccumData',deps_long,'Fun',@mean);

r=sqrt(dR_bins{1});


eulerStats.sigmaV = sigmaV;
eulerStats.r = r;
eulerStats.Nsamples = Nsamples;

eulerStats.S2x = S2x;
eulerStats.S2y = S2y;
eulerStats.S2z = S2z;

eulerStats.Splong{1} = S1long;
eulerStats.Splong{2} = S2long;
eulerStats.Splong{3} = S3long;
eulerStats.Splong{4} = S4long;
eulerStats.Splong{5} = S5long;

eulerStats.SplongAbs{1} = S1longAbs;
eulerStats.SplongAbs{2} = S2long;
eulerStats.SplongAbs{3} = S3longAbs;
eulerStats.SplongAbs{4} = S4long;
eulerStats.SplongAbs{5} = S5longAbs;

eulerStats.Sau = Sau;
eulerStats.Saulong = Sau_long;


%%
%figure;
%loglog(r,S2x*F_Ech^2*px^2);
%loglog(r,S2y*F_Ech^2*px^2);
%semilogx(r,Sau_long);
%hold on;
%loglog(r,S3long*F_Ech^3*px^3);
%loglog(r,S4long*F_Ech^4*px^4);
%loglog(r,S5long*F_Ech^5*px^5);
%plot(r,S2x*30^2*px^2);
%plot(r,S2y*30^2*px^2);
%semilogx(r,S3long);
%loglog(r,S4);
%loglog(r,S5);

%plot(r,r.^(2/3)*30,'--');

%plot(r,r.^(3/3)*450,'--');
%plot(r,r.^(4/3)*8000,'--');
%plot(r,r.^(5/3)*160000,'--');


%% Eulerian correlation function
NbinsSp = 2*Nbins;
dRmaxSp = dRmax;

sigmaV=sqrt(mean(S2long(end-4:end)));
eulerStats.sigmaVb = sigmaV;

[S2longb dR_edges dR_bins]=histcn([sqrt(dR2')],linspace(-dRmaxSp/(NbinsSp-1)/2,dRmaxSp,NbinsSp),'AccumData',[abs(dV_long).^2],'Fun',@mean);
r=dR_bins{1};


Ruu=(1-S2longb/sigmaV.^2);

%Ruu=Ruu-mean(Ruu(end-20:end));
%Ruu=Ruu./Ruu(1);
Ruu_sym=[flipud(Ruu(1:end)) ; Ruu];
r_sym=[-flipud(r(1:end)') ; r'];

eulerStats.Ruur = r_sym;
eulerStats.Ruu = Ruu_sym;
%% Eulerian Spectrum
PSD = (fft(Ruu_sym));
PSD = PSD(1:numel(Ruu));
k=[0:numel(Ruu)-1]*2*pi./max(r);
kresampled=linspace(0,max(k),512);
PSDresampled=interp1(k,abs(PSD),kresampled);

eulerStats.PSDk = kresampled;
eulerStats.PSD = PSDresampled;
% %%
% figure;loglog(kresampled,abs(PSDresampled))
% hold on
% loglog(kresampled,kresampled.^(-5/3),'--')
% axis([1e-2 2 1e-2 1e2])
% grid minor


