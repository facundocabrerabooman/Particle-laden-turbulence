function [Sn,L]=structFunc(sig,n,varargin)

if nargin>2
	dt=varargin{1};
else
	dt=(1:numel(sig))-1;
end

inc=increment(sig,dt);
inc=cell2struct(inc,'inc',2);

L=arrayfun(@(x)(numel(x.inc)),inc);
Sn=arrayfun(@(x)(nanmean(x.inc.^n)),inc);

end

function inc=increment(sig,varargin)

if nargin>1
	dt=varargin{1};
else
	dt=(1:numel(sig))-1;
end

inc=cell(numel(dt),1);
for k=1:numel(dt)
	ind=1:(numel(sig)-dt(k));
	inc{k}=sig(ind+dt(k))-sig(ind);
end

end