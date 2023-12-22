function [St]=Stokes_Mei(phi,gam,varargin)

% St=Stokes_Mei(phi,gam,varagin)
%
% calculates Stokes number based on the definition in:
%
% R. Mei, Velocity fidelity of flow racer particles
% Experiments in Fluids 22,pp. 1-13 (1996)
% 
% 1+3/16Re					Re<0.01
% 1+0.1315*Re^(0.82-0.05w)	0.01<Re<20
% 1+0.1935*Re^0.6302			20<Re<260
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

Re

% if Re<=0.01
% 	F=1+3/16*Re;
% elseif and((Re > 0.01),(Re <=20))
% 	n=0.82-0.05*log10(Re);
% 	F=1+0.1315*Re.^n;
% elseif and((Re>20), (Re <=260))
% 	F=1+0.1935*Re.^0.6302;
% else
% 	F=0;
% end
%n=0.82-0.2/3*log10(phi);

if Re<=2e5
    F=(1+0.150*Re.^0.681)+Re/24*0.407./(1+8710./Re);
else
    F=0;
end

St=1/36*phi.^2.*(1+2*gam)./F;

flag

if flag==1
    I=find(phi>1);
    St(I)=St(I).*phi.^(-2/3);
end