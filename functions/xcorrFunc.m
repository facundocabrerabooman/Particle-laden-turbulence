function [Sn ll s sdt]=xcorrFunc(sig,n,varargin)
if nargin>2
	dt=varargin{1};
else
	dt=(1:numel(sig))-1;
end

[inc, s, sdt] =increment(sig,dt,n);
inc=cell2struct(inc,'inc');


ll=arrayfun(@(x)(numel(x.inc)),inc);
%ll=ones(size(inc));
Sn=arrayfun(@(x)(mean(x.inc.^(n))),inc);


function [incr s sdt]=increment(sig,varargin)

%[incr]=increments_velocities(vel_in)
%
% signal des increments de vel_in : incr=vel_in(t+dt)=vel_in(t)

if nargin>2
    n = varargin{2};
    dt=varargin{1};
elseif nargin >1
    n=1;
	dt=varargin{1};
else
    n=1;
	dt=(1:numel(sig))-1;
end
k=1;
for kk=1:numel(dt)
	ind=1:(numel(sig)-dt(kk));
	incr{k}=(sig(ind+dt(kk)).*sig(ind));
    sdt(k) = mean(sig(ind+dt(kk)).^(2*n));
    s(k) = mean(sig(ind).^(2*n));
    k=k+1;
end



