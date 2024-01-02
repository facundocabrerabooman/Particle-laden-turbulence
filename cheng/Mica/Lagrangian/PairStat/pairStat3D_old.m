function pair = pairStat3D(tracks)
%% rearrange data "per image" instead of "per track"
%part=track2part(vtracks);

T=zeros(1,numel(tracks.xt));
Ntrack=zeros(1,numel(tracks.xt));

for k=1:numel(tracks.nstart)
        T(tracks.nstart(k):tracks.nend(k))=[tracks.tstart(k):tracks.tstop(k)];
        Ntracks(tracks.nstart(k):tracks.nend(k))=k;
end


%% Compute pair distances for all images
for k=1:numel(Tframe)
    I=find(T==Tframe(k));
    NN=numel(I);
    Ntrack=Ntracks(I);
    X=tracks.xt(I);
    Y=tracks.yt(I);
    Z=tracks.zt(I);
    Vx=tracks.dxt(I);
    Vy=tracks.dyt(I);
    Vz=tracks.dzt(I);
    if isfield(tracks,'dxpt');
        Vpx=tracks.dxpt(I);
        Vpy=tracks.dypt(I);
        Vpz=tracks.dzpt(I);
    end
    if isfield(tracks,'d2xt');
        Ax=tracks.d2xt(I);
        Ay=tracks.d2yt(I);
        Az=tracks.d2zt(I);
    end
    
    
    
    % Compute relative separations
    if isrow(X)
        X1=ones(NN,1)*X;
        X2=X'*ones(1,NN);
    else
        X1=X*ones(1,NN);
        X2=ones(NN,1)*X';
    end
    %dX2=(X2-X1).^2;
    dX=reshape(triu((X2-X1)),1,[]);
    Xm=reshape(triu((X2+X1)/2),1,[]);
   
    if isrow(Y)
        Y1=ones(NN,1)*Y;
        Y2=Y'*ones(1,NN);
    else
        Y1=Y*ones(1,NN);
        Y2=ones(NN,1)*Y';
    end
    %dY2=(Y2-Y1).^2;
    dY=reshape(triu((Y2-Y1)),1,[]);
    Ym=reshape(triu((Y2+Y1)/2),1,[]);
    
    if isrow(Z)
        Z1=ones(NN,1)*Z;
        Z2=Z'*ones(1,NN);
    else
        Z1=Z*ones(1,NN);
        Z2=ones(NN,1)*Z';
    end
    %dY2=(Y2-Y1).^2;
    dZ=reshape(triu((Z2-Z1)),1,[]);
    Zm=reshape(triu((Z2+Z1)/2),1,[]);
    
    dD2=(dX.^2+dY.^2+dZ.^2);
    
    % Compute relative velocities
    if isrow(Vx)
        X1=ones(NN,1)*Vx;
        X2=Vx'*ones(1,NN);
    else
        X1=Vx*ones(1,NN);
        X2=ones(NN,1)*Vx';
    end
    %dX2=(X2-X1).^2;
    dVx=reshape(triu((X2-X1)),1,[]);
    Vxm=reshape(triu((X2+X1)/2),1,[]);
    
    if isrow(Vy)
        Y1=ones(NN,1)*Vy;
        Y2=Vy'*ones(1,NN);
    else
        Y1=Vy*ones(1,NN);
        Y2=ones(NN,1)*Vy';
    end
    %dY2=(Y2-Y1).^2;
    dVy=reshape(triu((Y2-Y1)),1,[]);
    Vym=reshape(triu((Y2+Y1)/2),1,[]);
    
    if isrow(Vz)
        Z1=ones(NN,1)*Vz;
        Z2=Vz'*ones(1,NN);
    else
        Z1=Vz*ones(1,NN);
        Z2=ones(NN,1)*Vz';
    end
    %dY2=(Y2-Y1).^2;
    dVz=reshape(triu((Z2-Z1)),1,[]);
    Vzm=reshape(triu((Z2+Z1)/2),1,[]);
    
    dV2=(dVx.^2+dVy.^2+dVz.^2);
    V2m=(Vxm.^2+Vym.^2+Vzm.^2);
    
    % Compute relative fluctuating velocities
    if isfield(tracks,'dxpt');
        if isrow(Vpx)
            X1=ones(NN,1)*Vpx;
            X2=Vpx'*ones(1,NN);
        else
            X1=Vpx*ones(1,NN);
            X2=ones(NN,1)*Vpx';
        end
        %dX2=(X2-X1).^2;
        dVpx=reshape(triu((X2-X1)),1,[]);
        Vpxm=reshape(triu((X2+X1)/2),1,[]);
        
        if isrow(Vpy)
            Y1=ones(NN,1)*Vpy;
            Y2=Vpy'*ones(1,NN);
        else
            Y1=Vpy*ones(1,NN);
            Y2=ones(NN,1)*Vpy';
        end
        %dY2=(Y2-Y1).^2;
        dVpy=reshape(triu((Y2-Y1)),1,[]);
        Vpym=reshape(triu((Y2+Y1)/2),1,[]);
        
        if isrow(Vz)
            Z1=ones(NN,1)*Vpz;
            Z2=Vpz'*ones(1,NN);
        else
            Z1=Vpz*ones(1,NN);
            Z2=ones(NN,1)*Vpz';
        end
        %dY2=(Y2-Y1).^2;
        dVpz=reshape(triu((Z2-Z1)),1,[]);
        Vpzm=reshape(triu((Z2+Z1)/2),1,[]);
        
        dVp2=(dVpx.^2+dVpy.^2+dVpz.^2);
        Vp2m=(Vpxm.^2+Vpym.^2+Vpzm.^2);
    end
    
    % Compute relative accelerations
    if isfield(tracks,'d2xt');
        if isrow(Ax)
            X1=ones(NN,1)*Ax;
            X2=Ax'*ones(1,NN);
        else
            X1=Ax*ones(1,NN);
            X2=ones(NN,1)*Ax';
        end
        %dX2=(X2-X1).^2;
        dAx=reshape(triu((X2-X1)),1,[]);
        Axm=reshape(triu((X2+X1)/2),1,[]);
        
        if isrow(Ay)
            Y1=ones(NN,1)*Ay;
            Y2=Ay'*ones(1,NN);
        else
            Y1=Ay*ones(1,NN);
            Y2=ones(NN,1)*Ay';
        end
        %dY2=(Y2-Y1).^2;
        dAy=reshape(triu((Y2-Y1)),1,[]);
        Aym=reshape(triu((Y2+Y1)/2),1,[]);
        
        if isrow(Az)
            Z1=ones(NN,1)*Az;
            Z2=Az'*ones(1,NN);
        else
            Z1=Az*ones(1,NN);
            Z2=ones(NN,1)*Az';
        end
        %dY2=(Y2-Y1).^2;
        dAz=reshape(triu((Y2-Y1)),1,[]);
        Azm=reshape(triu((Z2+Z1)/2),1,[]);
        
        dA2=(dAx.^2+dAy.^2+dAz.^2);
        A2m=(Axm.^2+Aym.^2+Azm.^2);
    end
    %
    
    dNtrack1=ones(NN,1)*Ntrack;
    dNtrack2=Ntrack'*ones(1,NN);
    %dNtrack=dNtrack2+i*dNtrack1;
    dNtrack=reshape(triu(dNtrack2+i*dNtrack1),1,[]);
    
    II=find(dD2==0);
    dX(II)=[];
    dY(II)=[];
    dZ(II)=[];
    dD2(II)=[];
    Xm(II)=[];
    Ym(II)=[];
    Zm(II)=[];
    dVx(II)=[];
    dVy(II)=[];
    dVz(II)=[];
    dV2(II)=[];
    Vxm(II)=[];
    Vym(II)=[];
    Vzm(II)=[];
    V2m(II)=[];
    if isfield(tracks,'dxpt');
        dVpx(II)=[];
        dVpy(II)=[];
        dVpz(II)=[];
        dVp2(II)=[];
        Vpxm(II)=[];
        Vpym(II)=[];
        Vpzm(II)=[];
        Vp2m(II)=[];
    end
    if isfield(tracks,'d2xt');
        dAx(II)=[];
        dAy(II)=[];
        dAz(II)=[];
        dA2(II)=[];
        Axm(II)=[];
        Aym(II)=[];
        Azm(II)=[];
        A2m(II)=[];
    end
    dNtrack(II)=[];
    
    pair(k).dX=dX;
    pair(k).dY=dY;
    pair(k).dZ=dZ;
    pair(k).Xm=Xm;
    pair(k).Ym=Ym;
    pair(k).Zm=Zm;
    pair(k).dD2=dD2;
    pair(k).dVx=dVx;
    pair(k).dVy=dVy;
    pair(k).dVz=dVz;
    pair(k).dV2=dV2;
    pair(k).Vxm=Vxm;
    pair(k).Vym=Vym;
    pair(k).Vzm=Vzm;
    pair(k).V2m=V2m;
    if isfield(tracks,'dxpt');
        pair(k).dVpx=dVpx;
        pair(k).dVpy=dVpy;
        pair(k).dVpz=dVpz;
        pair(k).dVp2=dVp2;
        pair(k).Vpxm=Vpxm;
        pair(k).Vpym=Vpym;
        pair(k).Vpzm=Vpzm;
        pair(k).Vp2m=Vp2m;
    end
    if isfield(tracks,'d2xt');
        pair(k).dAx=dAx;
        pair(k).dAy=dAy;
        pair(k).dAz=dAz;
        pair(k).dA2=dA2;
        pair(k).Axm=Axm;
        pair(k).Aym=Aym;
        pair(k).Azm=Azm;
        pair(k).A2m=A2m;
    end
    pair(k).dNtrack=dNtrack;
    
end
 