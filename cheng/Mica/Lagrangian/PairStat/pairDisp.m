function meanD2 = pairDisp(vtracks,Dmax,Nbins)
%% rearrange data "per image" instead of "per track"
part=track2part(vtracks);
 
%% Compute pair distances for all images
for k=1:numel(part)
    NN=numel(part(k).X);
    
    % Compute relative separations
    if isrow(part(k).X)
        X1=ones(NN,1)*part(k).X;
        X2=part(k).X'*ones(1,NN);
    else
        X1=part(k).X*ones(1,NN);
        X2=ones(NN,1)*part(k).X';
    end
    %dX2=(X2-X1).^2;
    dX=reshape(triu((X2-X1)),1,[]);
   
    if isrow(part(k).Y)
        Y1=ones(NN,1)*part(k).Y;
        Y2=part(k).Y'*ones(1,NN);
    else
        Y1=part(k).Y*ones(1,NN);
        Y2=ones(NN,1)*part(k).Y';
    end
    %dY2=(Y2-Y1).^2;
    dY=reshape(triu((Y2-Y1)),1,[]);
    
    dR2=(dX.^2+dY.^2);
    
    % Compute relative velocities
    if isrow(part(k).Vx)
        X1=ones(NN,1)*part(k).Vx;
        X2=part(k).Vx'*ones(1,NN);
    else
        X1=part(k).Vx*ones(1,NN);
        X2=ones(NN,1)*part(k).Vx';
    end
    %dX2=(X2-X1).^2;
    dVx=reshape(triu((X2-X1)),1,[]);
   
    if isrow(part(k).Vy)
        Y1=ones(NN,1)*part(k).Vy;
        Y2=part(k).Vy'*ones(1,NN);
    else
        Y1=part(k).Vy*ones(1,NN);
        Y2=ones(NN,1) *part(k).Vy';
    end
    %dY2=(Y2-Y1).^2;
    dVy=reshape(triu((Y2-Y1)),1,[]);
    
    dV2=(dX.^2+dY.^2);

    % Compute relative accelerations
    if isrow(part(k).Ax)
        X1=ones(NN,1)*part(k).Ax;
        X2=part(k).Ax'*ones(1,NN);
    else
        X1=part(k).Ax*ones(1,NN);
        X2=ones(NN,1)*part(k).Ax';
    end
    %dX2=(X2-X1).^2;
    dAx=reshape(triu((X2-X1)),1,[]);
   
    if isrow(part(k).Ay)
        Y1=ones(NN,1)*part(k).Ay;
        Y2=part(k).Ay'*ones(1,NN);
    else
        Y1=part(k).Ay*ones(1,NN);
        Y2=ones(NN,1)*part(k).Ay';
    end
    %dY2=(Y2-Y1).^2;
    dAy=reshape(triu((Y2-Y1)),1,[]);
    
    dA2=(dX.^2+dY.^2);
    %
    
    dNtrack1=ones(NN,1)*part(k).Ntrack;
    dNtrack2=part(k).Ntrack'*ones(1,NN);
    %dNtrack=dNtrack2+i*dNtrack1;
    dNtrack=reshape(triu(dNtrack2+i*dNtrack1),1,[]);
    
    II=find(dD2==0);
    dX(II)=[];
    dY(II)=[];
    dR2(II)=[];
    dVx(II)=[];
    dVy(II)=[];
    dV2(II)=[];
    dAx(II)=[];
    dAy(II)=[];
    dAz(II)=[];
    dNtrack(II)=[];
    
    part(k).dX=dX;
    part(k).dY=dY;
    part(k).dR2=dR2;
    part(k).dVx=dVx;
    part(k).dVy=dVy;
    part(k).dV2=dV2;
    part(k).dAx=dAx;
    part(k).dAy=dAy;
    part(k).dA2=dA2;
    part(k).dNtrack=dNtrack;
    
    clear dX2 dY2 dD2 dNtrack II NN;
end
 
%% part2track for pairs
dNtrack=[part.dNtrack];
dNtrackU=unique(dNtrack);
dX=[part.dX];
dY=[part.dY];
dR2=[part.dR2];
dVx=[part.dVx];
dVy=[part.dVy];
dV2=[part.dV2];
dAx=[part.dAx];
dAy=[part.dAy];
dA2=[part.dA2];

 for k=1:numel(dNtrackU)
   II=find(dNtrack==dNtrackU(k));  
   trackpair(k).dX=dX(II);
   trackpair(k).dY=dY(II);
   trackpair(k).dR2=dR2(II);
   trackpair(k).dNtrack=dNtrackU(k);
end
 
%% define initial separations and binnings
%Nbins=64;
%Dmax=10^2; %you have to set this value to be the maximum separation you want to explore, now it is set for 100mm (assuming your data is in mm)
 
dX0=arrayfun(@(X)(abs(X.dX(1))),trackpair);
dY0=arrayfun(@(X)(abs(X.dY(1))),trackpair);
dR20=arrayfun(@(X)(X.dR2(1)),trackpair);
 
edgesX=linspace(0,Dmax,Nbins);
edgesY=linspace(0,Dmax,Nbins);
edgesR2=linspace(0,Dmax^2,Nbins);
 
binX=(edgesX(2:end)+edgesX(1:end-1))/2;
binY=(edgesY(2:end)+edgesY(1:end-1))/2;
binR2=(edgesR2(2:end)+edgesR2(1:end-1))/2;
 
%%
[NX,edgesX0,binX0]=histcounts(dX0,edgesX);
[NY,edgesY0,binY0]=histcounts(dY0,edgesY);
[NR2,edgesR20,binR20]=histcounts(dR20,edgesR2);
 
for k=2:numel(binX)    
 %   IX=find(binX0==k);
%    IY=find(binY0==k);
    IR=find(binR20==k); 
% for statistics on 1D separation, we still condition on the 2D separation 
% R20 to avoid combining for instance particles with small separation in X but actually very much separated in Y
    IX=IR; 
    IY=IR;
 
    dX2=arrayfun(@(X)((X.dX-X.dX(1)).^2),trackpair(IX),'UniformOutput',false);
    dY2=arrayfun(@(X)((X.dY-X.dY(1)).^2),trackpair(IY),'UniformOutput',false);
    dR2=arrayfun(@(X)(X.dR2-X.dR2(1)),trackpair(IR),'UniformOutput',false);   
    dD2=arrayfun(@(X)((X.dX-X.dX(1)).^2+(X.dY-X.dY(1)).^2),trackpair(IR),'UniformOutput',false);
 
    NX=cellfun(@(X)(numel(X)),dX2);
    NY=cellfun(@(X)(numel(X)),dY2);
    NR=cellfun(@(X)(numel(X)),dR2);
    ND=cellfun(@(X)(numel(X)),dD2);
 
    MX2=NaN(numel(dX2),max(NX));
    MY2=NaN(numel(dY2),max(NY));
    MR2=NaN(numel(dD2),max(NR));
    MD2=NaN(numel(dD2),max(ND));
 
 
    for kk=1:numel(dX2)
        MX2(kk,1:NX(kk))=dX2{kk};
    end
 
    for kk=1:numel(dY2)
        MY2(kk,1:NY(kk))=dY2{kk};
    end
 
    for kk=1:numel(dR2)
        MR2(kk,1:NR(kk))=dR2{kk};
    end
 
    for kk=1:numel(dD2)
        MD2(kk,1:ND(kk))=dD2{kk};
    end
 
    dX2=nanmean(MX2,1);
    dY2=nanmean(MY2,1);
    dR2=nanmean(MY2,1);
    dD2=nanmean(MD2,1);
 
    dX2(dX2==NaN)=[];
    dY2(dY2==NaN)=[];
    dR2(dR2==NaN)=[];
    dD2(dD2==NaN)=[];
    
    meanD2(k).dX2=nanmean(MX2,1);
    meanD2(k).dY2=nanmean(MY2,1);
    meanD2(k).dR2=nanmean(MR2,1);
    meanD2(k).dD2=nanmean(MD2,1);
    meanD2(k).TX=[0:numel(dX2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).TY=[0:numel(dY2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).TR=[0:numel(dR2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).TD=[0:numel(dD2)-1]/Fech;% you need to define Fech (sampling frequencey) in the workspace
    meanD2(k).dX20=mean(arrayfun(@(X)(X.dX(1).^2),trackpair(IR)));
    meanD2(k).dY20=mean(arrayfun(@(X)(X.dY(1).^2),trackpair(IR)));
    meanD2(k).dR20=binR2(k);
    clear MX2;
end
