function [df,idx]=finiteCenteredDiffn(f,dt,varargin);

% [df]=finiteDiffn(f,dt,varargin);
% df(:,1) contains the indexes for which finite differences could be
% calcultated
%
% calculates the centered of order o (o = varargin{1}, default o=1) of f.^p
% (p = varargin{2})
% sum from 0 to o  of (-1)^k pCk f(x+(p/2-i)*dt)
% for centered differences dt must be even.
%

o = 1;
p = 1;

if max(mod(dt,2)) ~=0
    disp('dt must be even');
    return
end

if nargin > 2
    o = varargin{1};
end
if nargin > 3
    p = varargin{2};
end

f = f.^p;

%% Relevant indexes common to all dt (Remove first and last indexes) 
N = max(size(f));
idx = [1 + o*max(dt)/2 : N - o*max(dt)/2];

df = zeros(numel(idx),numel(dt));

%% Calculates centered increments for all dt

for k = 1:numel(dt)
    for kk = 0:o
        df(:,k) = df(:,k)+(-1).^kk.*nchoosek(o,kk).*f(idx+(o/2-kk)*dt(k));
    end
end





