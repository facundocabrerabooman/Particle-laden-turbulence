
function [incr,idx]=increment(sig,varargin)

%[incr]=increments_velocities(vel_in)
%
% signal des increments de vel_in : incr=vel_in(t+dt)=vel_in(t)

if nargin>1
	dt=varargin{1};
else
	dt=(1:numel(sig))-1;
end

for kk=1:numel(dt)
	ind=1:numel(sig)-dt(kk);
    if min(ind)>0
        idx{kk}=ind+floor(dt(kk)/2);
        incr{kk}=(sig(ind+dt(kk))-sig(ind));
    end
end

