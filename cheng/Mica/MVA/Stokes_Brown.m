function [St]=Stokes_Brown(phi,gam,varargin)

% St=Stokes_Mei(phi,gam,varagin)
%
% calculates Stokes number based on the definition in:
%
% P. Brown and D. Lawler. Sphere drag and settling velocity revisited. Journal of Environmental Engineering-Asce, 129(3):222?231, Mar. 2003.
% 
% F=(1+0.150*Re.^0.681)+Re/24*0.407./(1+8710./Re);
%
% If varargin{1}==1 flow timescale is based on eddy turnover time at the particle
% dimension scale.
%
% varargin{2} can optionally specify a particle Reynolds number otherwise
% it is estimated as phi^(4/3) (i.e. based on fluid velocity at particle
% lengh scale).

flag=0;
if nargin>2
	flag=varargin{1};
end


Re=phi.^(4/3);

if nargin>3
    Re=varargin{2};
end

if Re<=2e5
    F=(1+0.150*Re.^0.681)+Re/24*0.407./(1+8710./Re);
else
    F=0;
end

St=1/36*phi.^2.*(1+2*gam)./F;
%St=1/18*phi.^2.*(gam)./F;
if flag==1
    I=find(phi>1);
    St(I)=St(I).*phi.^(-2/3);
end