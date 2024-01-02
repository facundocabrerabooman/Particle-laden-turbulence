function Y=Extrap_n_convs(X,kernel)
% Y=Extrap_n_convs(X,kernel) deactivate NaN and interpolate them
extrapPoints=13;


X=X(:)';

warning off

LL=length(X);
LKernel=fix((length(kernel))/2);

innerRange=1:LL;
%% remove NaNs
innerRange=innerRange(~isnan(X(:)));

%% interpolate them and extrapolate
Xi=interp1(innerRange,X(~isnan(X)),(-LKernel+1):(LL+LKernel),'pchip');

%% refine extrapolation using a cubic fit for the last 12 data points
use=min(extrapPoints,length(innerRange));

X=X(~isnan(X));

p=polyfit(innerRange(1:use),X(1:use),1);
Xi(1:LKernel)=polyval(p,(-LKernel+1):0);

p=polyfit(innerRange((end-use):end),X((end-use):end),1);
Xi((end-LKernel+1):end)=polyval(p,(LL+1):(LL+LKernel));

%% convolute
Y=conv(Xi,kernel,'valid');


warning on