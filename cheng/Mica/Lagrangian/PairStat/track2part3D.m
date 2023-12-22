function part=track2part(vtracks,varargin)

%N=arrayfun(@(X)(numel(X.Tf)),vtracks);
%N=arrayfun(@(X)(numel(X.frame)),vtracks);
N=arrayfun(@(X)(numel(X.Tf)),vtracks);

if nargin > 1
    T0 = varargin{1};
else
    T0 = 0;
end

T=[vtracks.Tf] + T0;

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
%Zf=Zf(I);
%Z=[vtracks.z];
%Z=vertcat(vtracks.z);


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
%Vz=Vz(I);

if isrow(vtracks(1).Ax)
    Ax=[vtracks.Ax];
else
    Ax=vertcat(vtracks.Ax);
end
%Ax=Ax(I);

if isrow(vtracks(1).Ay)
    Ay=[vtracks.Ay];
else
    Ay=vertcat(vtracks.Ay);
end
%Ay=Ay(I);

if isrow(vtracks(1).Az)
    Az=[vtracks.Az];
else
    Az=vertcat(vtracks.Az);
end
%Az=Az(I);




Ntrack=X*0;
NN=1;
for k=1:numel(vtracks)
    
    Ntrack(NN:NN+N(k)-1)=vtracks(k).Ntrack(1);
    
    T(NN:NN+N(k)-1)=T(NN:NN+N(k)-1);
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