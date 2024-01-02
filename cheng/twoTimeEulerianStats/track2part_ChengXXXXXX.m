function part=track2part_Cheng(vtracks,varargin)

if nargin<2
    
    %N=arrayfun(@(X)(numel(X.Tf)),vtracks);
    %N=arrayfun(@(X)(numel(X.frame)),vtracks);
    N=arrayfun(@(X)(numel(X.Tf)),vtracks);
    if isrow(vtracks(1).Tf)
        T=[vtracks.Tf];
    else
        T=vertcat(vtracks.Tf);
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
    
    for k=1:numel(binT)
        I=find(T==binT(k));
        part(k).T=binT(k);
        %   part(k).X=X(I);
        %   part(k).Y=Y(I);
        part(k).X=X(I);
        part(k).Y=Y(I);
        part(k).Z=Z(I);
        %   part(k).Theta=Theta(I);
        %   part(k).Rho=Rho(I);
        part(k).Vx=Vx(I);
        part(k).Vy=Vy(I);
        part(k).Vz=Vz(I);
        part(k).Ax=Ax(I);
        part(k).Ay=Ay(I);
        part(k).Az=Az(I);
        part(k).Ntrack=Ntrack(I);
    end
else
    

    fields = varargin{1};
    Tfield = fields{varargin{2}};
    
    N=arrayfun(@(X)(numel(X.(Tfield))),vtracks);
    
    
    if isrow(vtracks(1).(Tfield))
        T=[vtracks.(Tfield)];
    else
        T=vertcat(vtracks.(Tfield));
    end
    
    binT=unique(T);
    
    NN=1;
    for k=1:numel(vtracks)
        Ntrack(NN:NN+N(k)-1)=k;
        NN=NN+N(k);
    end
    
    folder = [pwd '\temp_part'];
    mkdir(folder)
    kk = 0;
    kp = 0;
    thres = 10;
    for k=1:numel(binT)
        kp = kp+1;
        part(kp).T=binT(k);
        I=find(T==binT(k));
        part(kp).Ntrack=Ntrack(I);
        
        for kfield=1:numel(fields)
            if isrow(vtracks(1).(fields{kfield}))
                X=[vtracks.(fields{kfield})];
            else
                X=vertcat(vtracks.(fields{kfield}));
            end
            part(kp).(fields{kfield})=X(I);
        end
 
        if kp>thres || k == numel(binT)
            kk = kk+1;
            disp([num2str(kk) '/' num2str(numel(binT)/thres)])
            save([folder '\temp_part_' num2str(kk) '.mat'],'part')
            clear part
            kp = 0;
        end
    end
    
    part = [];
    for kkk = 1:kk
        tp = load([folder '\temp_part_' num2str(kkk) '.mat']);
        part = [part;tp.part];
        clear tp
    end
end