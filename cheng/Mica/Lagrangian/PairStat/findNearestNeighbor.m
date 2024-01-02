function [dD,dX,dY,dZ] = findNearestNeighbor(V1,V2)

X1=V1(:,1);
Y1=V1(:,2);
Z1=V1(:,3);
X2=V2(:,1);
Y2=V2(:,2);
Z2=V2(:,3);

%% Compute pair distances
    
    % Compute relative separations
        NN1=numel(X1);
        NN2=numel(X2);
        X1=X1*ones(1,NN2);
        X2=ones(NN1,1)*X2';
        Y1=Y1*ones(1,NN2);
        Y2=ones(NN1,1)*Y2';
        Z1=Z1*ones(1,NN2);
        Z2=ones(NN1,1)*Z2';
        dX=triu((X2-X1));
        dY=triu((Y2-Y1));
        dZ=triu((Z2-Z1));
    
    dD=sqrt(dX.^2+dY.^2+dZ.^2);
    
 