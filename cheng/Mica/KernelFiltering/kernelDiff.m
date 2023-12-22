function data=kernelDiff(data,L,w,order)
%  data=kernelDiff(data,L,w,order) test function to try smoothed
%  derivatives
% 0 does olny smoothing


  %% computing   
if order==2
  kernel=accfiltcoef(w,L);
elseif order==1
  kernel=velfiltcoef(w,L);
else
  kernel=posfiltcoef(w,L);
end

data=Extrap_n_convs(data,kernel);

end
%% functions
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
    
    function kernel = posfiltcoef(filterwidth, filterlen)
    % calculate the coefficients associated with the
      % position gaussiane filter
      %
      % input:
      %   filterwidth     --  filter width (w in Nicolas' paper)
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
    
   