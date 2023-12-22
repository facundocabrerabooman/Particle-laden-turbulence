function kernel = accfiltcoef(filterwidth, filterlen)
% calculate the coefficients associated with the
% uses the formulas from the Diplomarbeit based on sums
% acceleration (2nd order derivative) gaussian filter
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

Tau2 = sum( t.^2 .* (2*(t/w).^2-1).*exp(-(t/w).^2));
TNull = sum((2*(t/w).^2-1).*exp(-(t/w).^2));

B=2/(L*(L+1)*(2*L+1)/3 -(2*L+1)*Tau2/TNull);
A=-(2*L+1)*B/TNull;
kernel = A*(2*(t/w).^2-1).*exp(-(t/w).^2)+B;
end

%% old with integrals
%     function kernel = accfiltcoef(filterwidth, filterlen)
% % calculate the coefficients associated with the
% % acceleration (2nd order derivative) gaussian filter
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
% fA = @(x,L,w)(x.^2-L^2/3).*(2*(x/w).^2-1).*exp(-(x/w).^2);
% fB = @(x,L,w)(2*(x/w).^2-1).*exp(-(x/w).^2);
%
% A=2/quad(@(x)fA(x,L,w),-L,L);
% B=-A*quad(@(x)fB(x,L,w),-L,L)/2/L;
%
% kernel = A*(2*(t/w).^2-1).*exp(-(t/w).^2)+B;