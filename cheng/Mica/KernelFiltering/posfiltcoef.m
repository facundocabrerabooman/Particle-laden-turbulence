function kernel = posfiltcoef(filterwidth, filterlen)
% calculate the coefficients associated with the
% position gaussiane filter
%
% input:
%   filterwidth --  filter width (w in Nicolas' paper)
%   filterlen   -- filter length (2*L+1 in Nicolas' paper)
%
% outputs:
%   kernel   --  coefficients of the filter, the filtered position is then
%           x_filt = sum(kernel.*x)
%
N = filterlen;
L = floor(N/2);
w = filterwidth;
t = (-L:L)';
A=1/sum(exp(-(t/w).^2));
kernel = A*exp(-(t/w).^2);

end

%% old using integrals:

% function kernel = posfiltcoef(filterwidth, filterlen)
% % calculate the coefficients associated with the
% % position gaussiane filter
% %
% % input:
% %   filterwidth     --  filter width (w in Nicolas' paper)
% %   filterlen   -- filter length (2*L+1 in Nicolas' paper)
% %
% % outputs:
% %   kernel   --  coefficients of the filter, the filtered position is then
% %           x_filt = sum(kernel.*x)
% %
% N = filterlen;
% L = floor(N/2);
% w = filterwidth;
% t = [-L:L]';
% B=exp(-L^2/w^2)/(2*L*exp(-L^2/w^2)-w*sqrt(pi)*erf(L/w));
% A=(1-2*B*L)/(w*sqrt(pi)*erf(L/w));
% kernel = A*exp(-(t/w).^2)+B;

