F_Ech=12000;
%%
dNtrack=[pair.dNtrack];
dX=[pair.dX];
dY=[pair.dY];
dD2=[pair.dD2];
dVx=[pair.dVx];
dVy=[pair.dVy];
dV2=[pair.dV2];
if isfield(pair,'dVpx');
    dVpx=[pair.dVpx];
    dVpy=[pair.dVpy];
    dVp2=[pair.dVp2];
end
dAx=[pair.dAx];
dAy=[pair.dAy];
dA2=[pair.dA2];
%Xm=[pair.Xm];
%Ym=[pair.Ym];
%%

[dNs,Is]=sort(dNtrack);
pairTracks = bwconncomp(1-abs(diff(dNs))>=1);
I0=cellfun(@(X)(X(1)),pairTracks.PixelIdxList);
[dNtrackU]=unique(dNtrack(Is));
dX=dX(Is);
dY=dY(Is);
dD2=dD2(Is);
dVx=dVx(Is);
dVy=dVy(Is);
dV2=dV2(Is);
if isfield(pair,'dVpx');
    dVpx=dVpx(Is);
    dVpy=dVpy(Is);
    dVpz=dVpz(Is);
    dVp2=dVp2(Is);
end
dAx=dAx(Is);
dAy=dAy(Is);
dAz=dAz(Is);
dA2=dA2(Is);
%%
%Xm=Xm(Is);
%Ym=Ym(Is);
%%
%[Thetam,Rhom,Zm]=cart2pol((Xm-nanmean(Xm))*90e-6,(Ym-nanmean(Ym))*90e-6,(Zm-nanmean(Zm))*90e-6);

%%
Ltrack=cellfun(@numel,pairTracks.PixelIdxList);
Lmax=max(Ltrack);
MD2=NaN(pairTracks.NumObjects,Lmax);

for kk=1:pairTracks.NumObjects
    MD2(kk,1:Ltrack(kk)+1)=...
        (dX(I0(kk):I0(kk)+Ltrack(kk))-dX(I0(kk))).^2+...
        (dY(I0(kk):I0(kk)+Ltrack(kk))-dY(I0(kk))).^2;
end

%%
D2max=650^2;
V2max=10;
A2max=3;

Nbins=32;

edgesD2=linspace(0,D2max,Nbins);
edgesV2=logspace(-2,log10(V2max),Nbins);
%edgesV2=linspace(0,V2max,Nbins);
edgesA2=linspace(0,A2max,Nbins);
%edgesRhom=linspace(0,max(Rhom)*1.0001,Nbins);

binD2=(edgesD2(2:end)+edgesD2(1:end-1))/2;
binV2=(edgesV2(2:end)+edgesV2(1:end-1))/2;
binA2=(edgesA2(2:end)+edgesA2(1:end-1))/2;
%binRhom=(edgesRhom(2:end)+edgesRhom(1:end-1))/2;
%binZm=(edgesZm(2:end)+edgesZm(1:end-1))/2;
%%
[ND2,edgesD20,binD20]=histcounts(dD2(I0),edgesD2);
[NV2,edgesV20,binV20]=histcounts(dV2(I0),edgesV2);
[NA2,edgesA20,binA20]=histcounts(dA2(I0),edgesA2);
%[NRhom,edgesRhom0,binRhom0]=histcounts(Rhom(I0),edgesRhom);
%[NZm,edgesZm0,binZm0]=histcounts(Zm(I0),edgesZm);
%% D2 conditionned on initial separation
for k=1:numel(binD2)
    I=find(binD20==k);
    meanD2(k).dD2=nanmean(MD2(I,:),1);
    meanD2(k).dD20=binD2(k);
    meanD2(k).dT=([1:Lmax+1]-1)/12000;
    meanD2(k).ND20=numel(I);
end

%% D2 conditionned on initial velocity
for k=1:numel(binV2)
    I=find(binV20==k);
    meanD2(k).dV2=nanmean(MD2(I,:),1);
    meanD2(k).dV20=binV2(k);
    meanD2(k).NV20=numel(I);
end

%% D2 conditionned on initial acceleration
for k=1:numel(binA2)
    I=find(binA20==k);
    meanD2(k).dA2=nanmean(MD2(I,:),1);
    meanD2(k).dA20=binA2(k);
    meanD2(k).NA20=numel(I);
end

%%
%% D2 conditionned on initial separation and velocity
I=find((binD20==2)&(binV20==1));
figure;loglog(meanD2(1).dT,nanmean(MD2(I,:),1));

%% D2 conditionned on initial velocity and position
%I=find(((binRhom0<=3)&(abs(binZm0-16)<=4))&(binV20<=2));
%figure;loglog(meanD2(1).dT,nanmean(MD2(I,:),1)./meanD2(1).dT.^3);