function [R ll]=autocorr(sig,moy,varargin)
if nargin>2
	dt=varargin{1};
else
	dt=(1:numel(sig))-1;
end


r=corr(sig-moy,dt);
r=cell2struct(r,'R');

ll=arrayfun(@(x)(numel(x.R)),r);

R=arrayfun(@(x)(mean(x.R)),r)';

%if nargin==3
	R=R/R(1);
%end


function [acorr]=corr(sig,varargin)

%[incr]=increments_velocities(vel_in)
%
% signal des increments de vel_in : incr=vel_in(t+dt)=vel_in(t)

if nargin>1
	dt=varargin{1};
else
	dt=(1:numel(sig))-1;
end

for kk=1:numel(dt)
	ind=1:(numel(sig)-dt(kk));
	acorr{kk}=sig(ind+dt(kk)).*sig(ind);
end