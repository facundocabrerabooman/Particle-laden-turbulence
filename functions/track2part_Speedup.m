function part=track2part_Speedup(vtracks,varargin)
% FCB changed  Tf to Tf sept 19th


folder = pwd;
% disp(pwd)

if nargin<2
    
    %N=arrayfun(@(X)(numel(X. Tf)),vtracks);s
    %N=arrayfun(@(X)(numel(X.frame)),vtracks);
    N=arrayfun(@(X)(numel(X. Tf)),vtracks);
    %if isrow(vtracks(1). Tf)
    if isrow(vtracks(1). Tf) 
        T=[vtracks. Tf];
    else
        T=vertcat(vtracks. Tf);
    end
    %T=vertcat(vtracks.frame);
    %[T,I]=sort(T);
    
    %X=[vtracks.X];
    %X=X(I);
    
    %Y=[vtracks.Y];
    %Y=Y(I);
    if isrow(vtracks(1).Xf)
        X=[vtracks.Xf];
    else
        X=vertcat(vtracks.Xf);
    end
    %Xf=Xf(I);
    %X=[vtracks.x];
    %X=vertcat(vtracks.x);
    
    if isrow(vtracks(1).Yf)
        Y=[vtracks.Yf];
    else
        Y=vertcat(vtracks.Yf);
    end
    %Yf=Yf(I);
    %Y=[vtracks.y];
    %Y=vertcat(vtracks.y);

    if isrow(vtracks(1).Zf)
        Z=[vtracks.Zf];
    else
        Z=vertcat(vtracks.Zf);
    end
    
    %Theta=[vtracks.Yf];
    
    %Rho=[vtracks.Rhof];
    
    if isrow(vtracks(1).Vx)
        Vx=[vtracks.Vx];
    else
        Vx=vertcat(vtracks.Vx);
    end
    %Vx=Vx(I);
    
    if isrow(vtracks(1).Vy)
        Vy=[vtracks.Vy];
    else
        Vy=vertcat(vtracks.Vy);
    end
    %Vy=Vy(I);

    if isrow(vtracks(1).Vz)
        Vz=[vtracks.Vz];
    else
        Vz=vertcat(vtracks.Vz);
    end
    
    if isrow(vtracks(1).Ax)
        Ax=[vtracks.Ax];
    else
        Ax=vertcat(vtracks.Ax);
    end
    %Ay=Ay(I);
    
    if isrow(vtracks(1).Ay)
        Ay=[vtracks.Ay];
    else
        Ay=vertcat(vtracks.Ay);
    end
    %Ax=Ax(I);

    if isrow(vtracks(1).Az)
        Az=[vtracks.Az];
    else
        Az=vertcat(vtracks.Az);
    end
    
    Ntrack=X*0;
    
    NN=1;
    for k=1:numel(vtracks)
        Ntrack(NN:NN+N(k)-1)=k;
        NN=NN+N(k);
    end
    %Ntrack=Ntrack(I);
    binT=unique(T);
    
    mkdir([folder '\temp_part'])
    parfor k=1:numel(binT)
        save_temp_part2(folder,T,binT,k,X,Y,Z,Vx,Vy,Vz,Ax,Ay,Az,Ntrack)
    end

    disp('sum up ... temp_part')
    part = [];
    for k=1:numel(binT)
        load([folder '\temp_part\temp_part_' num2str(k) '.mat']);
        part(k)= tp2;
    end
    disp('deleting temp_part')
    rmdir([folder '\temp_part'],'s')


else
    fields = varargin{1};
     Tfield = fields{varargin{2}};
    
    N=arrayfun(@(X)(numel(X.( Tfield))),vtracks);
    
    
    if isrow(vtracks(1).( Tfield))
        T=[vtracks.( Tfield)];
    else
        T=vertcat(vtracks.( Tfield));
    end
    
    binT=unique(T);
    disp(numel(binT))
    
    NN=1;
    for k=1:numel(vtracks)
        Ntrack(NN:NN+N(k)-1)=k;
        NN=NN+N(k);
    end
    
    mkdir([folder '\temp_part'])
    parfor k=1:numel(binT)
        save_temp_part1(folder,vtracks,fields,T,Ntrack,binT,k)
    end
    
    disp('sum up ... temp_part')
    imax = 1e5;
    part(1:imax) = struct('T',nan,'Ntrack',nan, 'Tf',nan,'Xf',nan,'Yf',nan,'Zf',nan,'Vx',nan,'Vy',nan,'Vz',nan,'Ax',nan,'Ay',nan,'Az',nan);
    for k=1:numel(binT)
        load([folder '\temp_part\temp_part_' num2str(k) '.mat']);
        part(k)= tp1;
    end
    part = part(~isnan([part.T]));
    disp('deleting ... temp_part')
    rmdir([folder '\temp_part'],'s')
    save([folder '\EulerianPart.mat'],'part')
end
