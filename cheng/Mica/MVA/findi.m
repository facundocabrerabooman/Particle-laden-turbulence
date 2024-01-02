function [i]=findi(var,val)

% i1=min(find(var>=val));
% i2=max(find(var<=val));
% 
% if abs(var(i1)-val)<abs(var(i2)-val)
%   i=i1;
% else
%   i=i2;
% end
[i,vali,delta]=findi2(var,val,1);