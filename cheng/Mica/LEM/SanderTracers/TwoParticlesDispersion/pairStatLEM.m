function [N dR S2 S3 Sau] = pairStatLEM(tracks)

Nfiles = arrayfun(@(X)(abs(imag(X.Ntrack(1)))),tracks);

NfilesUnique = unique(Nfiles);

% I=find(Nfiles == 1);
% part = track2part3D(tracks(I));
% pair = pairStat3D(part);
% pair = addStructFun(pair,'dX','Nfile',@(X)(1));
% trackPair = part2track3D(pair);
% trackPair = addStructFun(trackPair,'dX','Nfile',@(X)(1));



drmin = .5; %mm
drmax = 90; %mm
nbins = 100;
dredges = logspace(log10(drmin),log10(drmax),nbins);

S2.x=zeros(nbins-1,1);
S2.y=zeros(nbins-1,1);
S2.z=zeros(nbins-1,1);
S2.long=zeros(nbins-1,1);
S3.long=zeros(nbins-1,1);
Sau.tot=zeros(nbins-1,1);
Sau.long=zeros(nbins-1,1);
N = zeros(nbins-1,1);

for kfile = 1:numel(NfilesUnique)
    kfile
    I=find(Nfiles == kfile);
    part = track2part3D(tracks(I));
    pair = pairStat3D(part);
    pair = addStructFun(pair,'dX','Nfile',@(X)(kfile));
    %trackPair = part2track3D(pair);
    %trackPair = addStructFun(trackPair,'dX','Nfile',@(X)(kfile));
    %pair=appendStruct(pair,pair_tmp);
    %    trackPair=appendStruct(trackPair,trackPair_tmp);
    %%
    dX=[pair.dX];
    dY=[pair.dY];
    dZ=[pair.dZ];
    dR2=[pair.dR2];
    
    dVx=[pair.dVx];
    dVy=[pair.dVy];
    dVz=[pair.dVz];
    dV2=[pair.dV2];
    
    dAx=[pair.dAx];
    dAy=[pair.dAy];
    dAz=[pair.dAz];
    dA2=[pair.dA2];
    
    %% Calculate transverse and longitudinal increments
    %[theta,r]=cart2pol(dX,dY);
    %dA_trans=-sin(theta).*dAx+cos(theta).*dAy;
    %dV_trans=-sin(theta).*dVx+cos(theta).*dVy;
    
    
    dA_long=(dAx.*dX + dAy.*dY + dAz.*dZ)./(sqrt(dR2));
    dV_long=(dVx.*dX + dVy.*dY + dVz.*dZ)./(sqrt(dR2));
    
    deps = dAx.*dVx + dAy.*dVy;
    deps_long=dA_long.*dV_long;
    
    %% Calculate structure functions
    
    %[S2long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',dV_long.^2,'Fun',@mean);
    %[S3x dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',dVx.^3,'Fun',@mean);
    %[S3y dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',dVy.^3,'Fun',@mean);
    [N_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges);
    [S2x_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges,'AccumData',abs(dVx).^2,'Fun',@sum);
    [S2y_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges,'AccumData',abs(dVy).^2,'Fun',@sum);
    [S2z_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges,'AccumData',abs(dVz).^2,'Fun',@sum);
    
    [S2long_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges,'AccumData',abs(dV_long).^2,'Fun',@sum);
    [S3long_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges,'AccumData',(dV_long).^3,'Fun',@sum);
    %    [S4long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',dV_long.^4,'Fun',@mean);
    %    [S5long dR_edges dR_bins]=histcn(sqrt(dR2'),logspace(log10(5),log10(250),40),'AccumData',abs(dV_long.^5),'Fun',@mean);
    
    [Sau_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges,'AccumData',deps,'Fun',@sum);
    [Saulong_tmp dR_edges dR_bins]=histcn(sqrt(dR2'),dredges,'AccumData',deps_long,'Fun',@sum);
    
    N = N + N_tmp;
    S2.x = S2.x + S2x_tmp;
    S2.y = S2.y + S2y_tmp;
    S2.z = S2.z + S2z_tmp;
    S2.long = S2.long + S2long_tmp;
    S3.long = S3.long + S3long_tmp;
    Sau.tot = Sau.tot + Sau_tmp;
    Sau.long = Sau.long + Saulong_tmp;
    
end

dR = dR_bins{1};

