function [pair,trackpair]=pairStat(part)
%part=pairStat(vtracks)
%
% For a given part structure, calculates relative separation,
% velocity and acceleration of particles frame by frame.
% If you have a track strucuture, use track2part first.



% part=track2part(vtracks);
doVel = 1;
doAcc = 1;
%%

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
    
    dR2=(dX.^2+dY.^2);
    
    %[cTh,cRho]=cart2pol(cX-700,cY-372); %% polar coordinates of COM of pairs
    
    if doVel == 1
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
        
        dV2=(dVx.^2+dVy.^2);
    end
    %% Compute relative accelerations
    if doAcc==1
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
        
        dA2=(dAx.^2+dAy.^2);
    end
    %%
    
    if isrow(part(k).Ntrack)
        dNtrack1=ones(NN,1)*part(k).Ntrack;
        dNtrack2=part(k).Ntrack'*ones(1,NN);
    else
        dNtrack1=part(k).Ntrack*ones(1,NN);
        dNtrack2=ones(NN,1)*part(k).Ntrack';
    end
    dNtrack=dNtrack2+i*dNtrack1;
    dNtrack=reshape(triu(dNtrack,1),1,[]);
    
    II=find(abs(dNtrack==0));
    dX(II)=[];
    dY(II)=[];
    dR2(II)=[];
    
    if doVel == 1
        dVx(II)=[];
        dVy(II)=[];
        dV2(II)=[];
    end
    
    if doAcc == 1
        dAx(II)=[];
        dAy(II)=[];
        dA2(II)=[];
    end
    
    dNtrack(II)=[];
    %cRho(II)=[];
    
    pair(k).dX=dX;
    pair(k).dY=dY;
    pair(k).dR2=dR2;
    
    if doVel == 1
        pair(k).dVx=dVx;
        pair(k).dVy=dVy;
        pair(k).dV2=dV2;
    end
    
    if doAcc == 1
        pair(k).dAx=dAx;
        pair(k).dAy=dAy;
        pair(k).dA2=dA2;
    end
    
    pair(k).dNtrack=dNtrack;
    % pair(k).cRho=cRho;
    k
end

% %% part2track for pairs
% dNtrack=[part.dNtrack];
% dNtrackU=unique(dNtrack);
% dX=[part.dX];
% dY=[part.dY];
% dR2=[part.dR2];
% dVx=[part.dVx];
% dVy=[part.dVy];
% dV2=[part.dV2];
% %dAx=[part.dAx];
% %dAy=[part.dAy];
% %dA2=[part.dA2];
%
% for k=1:numel(dNtrackU)
%    II=find(dNtrack==dNtrackU(k));
%    trackpair(k).dX=dX(II);
%    trackpair(k).dY=dY(II);
%    trackpair(k).dR2=dR2(II);
%    trackpair(k).dVx=dVx(II);
%    trackpair(k).dVy=dVy(II);
%    trackpair(k).dV2=dV2(II);
%    trackpair(k).dNtrack=dNtrackU(k);
% end

