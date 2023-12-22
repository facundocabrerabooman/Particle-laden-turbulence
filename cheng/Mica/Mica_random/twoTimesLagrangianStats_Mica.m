function lagStats = twoTimesLagrangianStats_Mica(vtracks,varargin)
if nargin > 1
    Nmax = varargin{1};
end

sigmaR2 = mean((reshape([vtracks.Rhof],1,[])).^2);
MSDx = structFunc_struct(vtracks,'X',2,Nmax);
MSDy = structFunc_struct(vtracks,'Y',2,Nmax);

lagStats.MSDx = MSDx;
lagStats.MSDy = MSDy;
lagStats.sigmaR2 = sigmaR2;
