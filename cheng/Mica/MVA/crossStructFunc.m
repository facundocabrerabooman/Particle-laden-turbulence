function [Sn ll]=crossStructFunc(sig1,sig2,n,varargin)
if nargin>3
	dt=varargin{1};
else
	dt=(1:numel(sig1))-1;
end

inc1=increment(sig1,dt);
inc1=cell2struct(inc1,'inc');

inc2=increment(sig2,dt);
inc2=cell2struct(inc2,'inc');


ll=arrayfun(@(x)(numel(x.inc)),inc1);

Sn=arrayfun(@(x,y)(mean((x.inc.*y.inc).^n)),inc1,inc2);


function [incr]=increment(sig,varargin)

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
	incr{kk}=sig(ind+dt(kk))-sig(ind);
end

