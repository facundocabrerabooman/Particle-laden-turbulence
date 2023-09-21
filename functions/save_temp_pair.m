function save_temp_pair(folder,part,k)
    
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
    dX=reshape(triu((X2-X1),1),1,[]);
    cX=reshape(triu((X2+X1),1)/2,1,[]);
   
    if isrow(part(k).Y)
        Y1=ones(NN,1)*part(k).Y;
        Y2=part(k).Y'*ones(1,NN);
    else
        Y1=part(k).Y*ones(1,NN);
        Y2=ones(NN,1)*part(k).Y';
    end
    %dY2=(Y2-Y1).^2;
    dY=reshape(triu((Y2-Y1),1),1,[]);
    cY=reshape(triu((Y2+Y1),1)/2,1,[]);
    
    if isrow(part(k).Z)
        Z1=ones(NN,1)*part(k).Z;
        Z2=part(k).Z'*ones(1,NN);
    else
        Z1=part(k).Z*ones(1,NN);
        Z2=ones(NN,1)*part(k).Z';
    end
    %dZ2=(Z2-Z1).^2;
    dZ=reshape(triu((Z2-Z1),1),1,[]);
    cZ=reshape(triu((Z2+Z1),1)/2,1,[]);
    
    dR2=(dX.^2+dY.^2+dZ.^2);
    
    %[cTh,cRho]=cart2pol(cX-700,cY-372); %% polar coordinates of COM of pairs

    %% Compute relative velocities
    if isrow(part(k).Vx)
        X1=ones(NN,1)*part(k).Vx;
        X2=part(k).Vx'*ones(1,NN);
    else
        X1=part(k).Vx*ones(1,NN);
        X2=ones(NN,1)*part(k).Vx';
    end
    %dX2=(X2-X1).^2;
    dVx=reshape(triu((X2-X1),1),1,[]);
   
    if isrow(part(k).Vy)
        Y1=ones(NN,1)*part(k).Vy;
        Y2=part(k).Vy'*ones(1,NN);
    else
        Y1=part(k).Vy*ones(1,NN);
        Y2=ones(NN,1) *part(k).Vy';
    end
    %dY2=(Y2-Y1).^2;
    dVy=reshape(triu((Y2-Y1),1),1,[]);
    
    if isrow(part(k).Vz)
        Z1=ones(NN,1)*part(k).Vz;
        Z2=part(k).Vz'*ones(1,NN);
    else
        Z1=part(k).Vz*ones(1,NN);
        Z2=ones(NN,1) *part(k).Vz';
    end
    %dZ2=(Z2-Z1).^2;
    dVz=reshape(triu((Z2-Z1),1),1,[]);
    
    dV2=(dVx.^2+dVy.^2+dVz.^2);

    %% Compute relative accelerations
    if isrow(part(k).Ax)
        X1=ones(NN,1)*part(k).Ax;
        X2=part(k).Ax'*ones(1,NN);
    else
        X1=part(k).Ax*ones(1,NN);
        X2=ones(NN,1)*part(k).Ax';
    end
    %dX2=(X2-X1).^2;
    dAx=reshape(triu((X2-X1),1),1,[]);
   
    if isrow(part(k).Ay)
        Y1=ones(NN,1)*part(k).Ay;
        Y2=part(k).Ay'*ones(1,NN);
    else
        Y1=part(k).Ay*ones(1,NN);
        Y2=ones(NN,1)*part(k).Ay';
    end
    %dY2=(Y2-Y1).^2;
    dAy=reshape(triu((Y2-Y1),1),1,[]);
    
    if isrow(part(k).Az)
        Z1=ones(NN,1)*part(k).Az;
        Z2=part(k).Az'*ones(1,NN);
    else
        Z1=part(k).Az*ones(1,NN);
        Z2=ones(NN,1)*part(k).Az';
    end
    %dY2=(Y2-Y1).^2;
    dAz=reshape(triu((Z2-Z1),1),1,[]);    
    
    dA2=(dAx.^2+dAy.^2+dAz.^2);
    
    %%
    if isrow(part(k).Ntrack)
        dNtrack1=ones(NN,1)*real(part(k).Ntrack);
        dNtrack2=real(part(k).Ntrack')*ones(1,NN);
    else
        dNtrack1=real(part(k).Ntrack)*ones(1,NN);
        dNtrack2=ones(NN,1)*real(part(k).Ntrack');
    end
    dNtrack=dNtrack2+i*dNtrack1;
    dNtrack=reshape(triu(dNtrack,1),1,[]);
    
    II=find(abs(dNtrack==0));
    dX(II)=[];
    dY(II)=[];
    dZ(II)=[];
    dR2(II)=[];
    dVx(II)=[];
    dVy(II)=[];
    dVz(II)=[];
    dV2(II)=[];
    dAx(II)=[];
    dAy(II)=[];
    dAz(II)=[];
    dA2(II)=[];
    dNtrack(II)=[];

    tp.dX=dX;
    tp.dY=dY;
    tp.dZ=dZ;
    tp.dR2=dR2;
    tp.dVx=dVx;
    tp.dVy=dVy;
    tp.dVz=dVz;
    tp.dV2=dV2;
    tp.dAx=dAx;
    tp.dAy=dAy;
    tp.dAz=dAz;
    tp.dA2=dA2;
    tp.dNtrack=dNtrack;

    save([folder filesep 'temp_pair' filesep 'temp_pair_' num2str(k) '.mat'],'tp')