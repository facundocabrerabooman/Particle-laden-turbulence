function [SnBwd SnFwd ll]=structFuncCentered(sig,n,varargin)
% [Sn ll]=SructFunc(S1,field1,n)
% 
% calculates Sn centered with forward and backwart time lag for a single
% track
%

k=1;
N = numel(sig);
iCenter = round(N/2);
sigFwd = sig(iCenter:end);
sigBwd = fliplr(sig(1:iCenter));

if nargin>2
	dt=varargin{1};
else
	dt=(1:min(numel(sigFwd),numel(sigBwd)))-1;
end

[incBwd,incFwd]=increment(sig,dt);
incBwd=cell2struct(incBwd,'inc');
incFwd=cell2struct(incFwd,'inc');


ll=arrayfun(@(x)(numel(x.inc)),incBwd);
%ll=ones(size(inc));
SnBwd=arrayfun(@(x)(mean(x.inc.^n)),incBwd);
SnFwd=arrayfun(@(x)(mean(x.inc.^n)),incFwd);

function [incrBwd,incrFwd]=increment(sig,varargin)

%[incr]=increments_velocities(vel_in)
%
% signal des increments de vel_in : incr=vel_in(t+dt)=vel_in(t)


k=1;
N = numel(sig);
iCenter = round(N/2);
sigFwd = sig(iCenter:end);
sigBwd = fliplr(sig(1:iCenter));

if nargin>1
	dt=varargin{1};
else
	dt=(1:min(numel(sigFwd),numel(sigBwd)))-1;
end

for kk=1:numel(dt)
	ind=1:(numel(sigFwd)-dt(kk));
	incrFwd{k}=sigFwd(ind+dt(kk))-sigFwd(ind);
    ind=1:(numel(sigBwd)-dt(kk));
    incrBwd{k}=sigBwd(ind+dt(kk))-sigBwd(ind);
    k=k+1;
end

