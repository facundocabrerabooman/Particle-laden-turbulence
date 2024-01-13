function kernel = velfiltcoef(filterwidth, filterlen)
% calculate the coefficients associated with the
% velocity (1st order derivative) gaussain filter
%
% input:
%   filterwidth     --  filter width (w in Nicolas' paper)
%   filterlen   -- filter length (2*L+1 in Nicolas' paper)
%
% outputs:
%   kernel   --  coefficients of the filter, the filtered velocity is then
%           velfilt = sum(kernel.*x)

N = filterlen;
L = floor(N/2);
w = filterwidth;
t = (-L:L)';

A= -1/sum(t.^2.*exp(-(t/w).^2));
kernel = A*t.*exp(-(t/w).^2);
end


% %% old using integrals
%
% function kernel = velfiltcoef(filterwidth, filterlen)
% % calculate the coefficients associated with the
% % velocity (1st order derivative) gaussain filter
% %
% % input:
% %   filterwidth     --  filter width (w in Nicolas' paper)
% %   filterlen   -- filter length (2*L+1 in Nicolas' paper)
% %
% % outputs:
% %   kernel   --  coefficients of the filter, the filtered velocity is then
% %           velfilt = sum(kernel.*x)
% N = filterlen;
% L = floor(N/2);
% w = filterwidth;
% t = [-L:L]';
%
% A = -2.0 / (w ^2 * (sqrt(pi) * w * erf(L / w) - 2.0 * L * exp(-L^2 / w^2)));
%
% kernel = A*t.*exp(-(t/w).^2);
%
