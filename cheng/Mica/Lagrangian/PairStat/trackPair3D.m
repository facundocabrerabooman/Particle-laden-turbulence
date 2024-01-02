function meanD2=trackPair3D(part,Dmax,Nbins,Fech)

%% part2track for pairs
dNtrack=[part.dNtrack];
dNtrackU=unique(dNtrack);
dX=[part.dX];
dY=[part.dY];
dZ=[part.dZ];
dR2=[part.dR2];
dVx=[part.dVx];
dVy=[part.dVy];
dVz=[part.dVz];
dV2=[part.dV2];
if isfield(part,'dVpx');
    dVpx=[part.dVpx];
    dVpy=[part.dVpy];
    dVpz=[part.dVpz];
    dVp2=[part.dVp2];
end
dAx=[part.dAx];
dAy=[part.dAy];
dAz=[part.dAz];
dA2=[part.dA2];

 h = waitbar(0,'Calculating pair statistics...',...
             'CreateCancelBtn',...
             'setappdata(gcbf,''canceling'',1)');
 setappdata(h,'canceling',0)

%h = waitbar(0,'Calculating pair statistics...');
Ntot=numel(dNtrackU)
for k=1:Ntot
   
   II=find(dNtrack==dNtrackU(k));  
   trackpair(k).dX=dX(II);
   trackpair(k).dY=dY(II);
   trackpair(k).dZ=dZ(II);
   trackpair(k).dR2=dR2(II);
   trackpair(k).dVx=dVx(II);
   trackpair(k).dVy=dVy(II);
   trackpair(k).dVz=dVz(II);
   trackpair(k).dV2=dV2(II);
   if isfield(part,'dVpx');
          trackpair(k).dVpx=dVpx(II);
          trackpair(k).dVpy=dVpy(II);
          trackpair(k).dVpz=dVpz(II);
          trackpair(k).dVp2=dVp2(II);
   end
   trackpair(k).dAx=dAx(II);
   trackpair(k).dAy=dAy(II);
   trackpair(k).dAz=dAz(II);
   trackpair(k).dA2=dA2(II);
   trackpair(k).dNtrack=dNtrackU(k);
 end
delete(h);
 
%% define initial separations and binnings
%Nbins=64;
%Dmax=10^2; %you have to set this value to be the maximum separation you want to explore, now it is set for 100mm (assuming your data is in mm)
 
dX0=arrayfun(@(X)(abs(X.dX(1))),trackpair);
dY0=arrayfun(@(X)(abs(X.dY(1))),trackpair);
dZ0=arrayfun(@(X)(abs(X.dZ(1))),trackpair);
dR20=arrayfun(@(X)(X.dR2(1)),trackpair);
 
edgesX=linspace(0,Dmax,Nbins);
edgesY=linspace(0,Dmax,Nbins);
edgesZ=linspace(0,Dmax,Nbins);
edgesR2=linspace(0,Dmax^2,Nbins);
 
binX=(edgesX(2:end)+edgesX(1:end-1))/2;
binY=(edgesY(2:end)+edgesY(1:end-1))/2;
binZ=(edgesY(2:end)+edgesY(1:end-1))/2;
binR2=(edgesR2(2:end)+edgesR2(1:end-1))/2;
 
%% stats conditionned to initial separation
[NX,edgesX0,binX0]=histcounts(dX0,edgesX);
[NY,edgesY0,binY0]=histcounts(dY0,edgesY);
[NY,edgesY0,binZ0]=histcounts(dZ0,edgesZ);
[NR2,edgesR20,binR20]=histcounts(dR20,edgesR2);
 
for k=2:numel(binX)    
%   IX=find(binX0==k);
%   IY=find(binY0==k);
%   IZ=find(binY0==k);
    IR=find(binR20==k); 
% for statistics on 1D separation, we still condition on the 2D separation 
% R20 to avoid combining for instance particles with small separation in X but actually very much separated in Y
    IX=IR; 
    IY=IR;
    IZ=IR;
 
    dX2=arrayfun(@(X)((X.dX-X.dX(1)).^2),trackpair(IX),'UniformOutput',false);
    dY2=arrayfun(@(X)((X.dY-X.dY(1)).^2),trackpair(IY),'UniformOutput',false);
    dZ2=arrayfun(@(X)((X.dZ-X.dZ(1)).^2),trackpair(IZ),'UniformOutput',false);
    dR2=arrayfun(@(X)(X.dR2-X.dR2(1)),trackpair(IR),'UniformOutput',false);   
    dD2=arrayfun(@(X)((X.dX-X.dX(1)).^2+(X.dY-X.dY(1)).^2),trackpair(IR),'UniformOutput',false);
 
    NX=cellfun(@(X)(numel(X)),dX2);
    NY=cellfun(@(X)(numel(X)),dY2);
    NZ=cellfun(@(X)(numel(X)),dZ2);
    NR=cellfun(@(X)(numel(X)),dR2);
    ND=cellfun(@(X)(numel(X)),dD2);
 
    MX2=NaN(numel(dX2),max(NX));
    MY2=NaN(numel(dY2),max(NY));
    MZ2=NaN(numel(dY2),max(NZ));
    MR2=NaN(numel(dD2),max(NR));
    MD2=NaN(numel(dD2),max(ND));
 
 
    for kk=1:numel(dX2)
        MX2(kk,1:NX(kk))=dX2{kk};
    end
 
    for kk=1:numel(dY2)
        MY2(kk,1:NY(kk))=dY2{kk};
    end
 
    for kk=1:numel(dZ2)
        MZ2(kk,1:NZ(kk))=dZ2{kk};
    end
    
    for kk=1:numel(dR2)
        MR2(kk,1:NR(kk))=dR2{kk};
    end
 
    for kk=1:numel(dD2)
        MD2(kk,1:ND(kk))=dD2{kk};
    end
 
    dX2=nanmean(MX2,1);
    dY2=nanmean(MY2,1);
    dZ2=nanmean(MZ2,1);
    dR2=nanmean(MY2,1);
    dD2=nanmean(MD2,1);
 
    dX2(isnan(dX2))=[];
    dY2(isnan(dY2))=[];
    dZ2(isnan(dZ2))=[];
    dR2(isnan(dR2))=[];
    dD2(isnan(dD2))=[];
    
    meanD2(k).dX2=nanmean(MX2,1);
    meanD2(k).dY2=nanmean(MY2,1);
    meanD2(k).dZ2=nanmean(MZ2,1);
    meanD2(k).dR2=nanmean(MR2,1);
    meanD2(k).dD2=nanmean(MD2,1);
    meanD2(k).TX=[0:numel(dX2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).TY=[0:numel(dY2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).TZ=[0:numel(dZ2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).TR=[0:numel(dR2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).TD=[0:numel(dD2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).dX20=mean(arrayfun(@(X)(X.dX(1).^2),trackpair(IR)));
    meanD2(k).dY20=mean(arrayfun(@(X)(X.dY(1).^2),trackpair(IR)));
    meanD2(k).dZ20=mean(arrayfun(@(X)(X.dZ(1).^2),trackpair(IR)));
    meanD2(k).dR20=binR2(k);
    clear MX2;
end
