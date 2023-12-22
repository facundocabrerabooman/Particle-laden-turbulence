function [neighborAll, neighborLayer] = neighborIdx(part,tracer,idx1,idx2,Rmin,Rmax,vargin)

dR = 0.5;

if nargin>6
    dR = vargin(1);
end


part0.xf = part.Xf(idx1);
part0.yf = part.Yf(idx1);
part0.zf = part.Zf(idx1);
part0.vx = part.Vx(idx1);
part0.vy = part.Vy(idx1);
part0.vz = part.Vz(idx1);
part0.ax = part.Ax(idx1);
part0.ay = part.Ay(idx1);
part0.az = part.Az(idx1);

tracer0.xf = tracer.Xf(idx2);
tracer0.yf = tracer.Yf(idx2);
tracer0.zf = tracer.Zf(idx2);
tracer0.vx = tracer.Vx(idx2);
tracer0.vy = tracer.Vy(idx2);
tracer0.vz = tracer.Vz(idx2);
tracer0.ax = tracer.Ax(idx2);
tracer0.ay = tracer.Ay(idx2);
tracer0.az = tracer.Az(idx2);

%% neighboring tracers

neighborAll = struct('idx',[],'d',[]);
neighborAll.Rmin = Rmin;
neighborAll.Rmax = Rmax;


for i = 1:numel(tracer0.xf)
    d(i,:) = sqrt((part0.xf-tracer0.xf(i))^2+(part0.yf-tracer0.yf(i))^2+(part0.zf-tracer0.zf(i))^2);
end

%% global neighboring
idx01 = find(d>Rmin & d<Rmax);
neighborAll.idx = idx2(idx01);
neighborAll.d = d(idx01);


%% front and back 
rpt.x = tracer0.xf(idx01) - part0.xf;
rpt.y = tracer0.yf(idx01) - part0.yf;
rpt.z = tracer0.zf(idx01) - part0.zf;

rv = part0.vx*rpt.x + part0.vy*rpt.y + part0.vz*rpt.z;
neighborAll.idxfront = neighborAll.idx(rv>0);
neighborAll.idxback = neighborAll.idx(rv<0);

%% global infos
% neighborAll.VrelFrontx = part0.vx - mean(tracer0.vx(idx01(rv>0)));
% neighborAll.VrelFronty = part0.vy - mean(tracer0.vy(idx01(rv>0)));
% neighborAll.VrelFrontz = part0.vz - mean(tracer0.vz(idx01(rv>0)));
% neighborAll.VrelFront = sqrt( (neighborAll.VrelFrontx)^2 + (neighborAll.VrelFronty)^2 + (neighborAll.VrelFrontz)^2);
% 
% neighborAll.VrelBackx = part0.vx - mean(tracer0.vx(idx01(rv<0)));
% neighborAll.VrelBacky = part0.vy - mean(tracer0.vy(idx01(rv<0)));
% neighborAll.VrelBackz = part0.vz - mean(tracer0.vz(idx01(rv<0)));
% neighborAll.VrelBack = sqrt( (neighborAll.VrelBackx)^2 + (neighborAll.VrelBacky)^2 + (neighborAll.VrelBackz)^2);
% 
% neighborAll.ArelFrontx = part0.ax - mean(tracer0.ax(idx01(rv>0)));
% neighborAll.ArelFronty = part0.ay - mean(tracer0.ay(idx01(rv>0)));
% neighborAll.ArelFrontz = part0.az - mean(tracer0.az(idx01(rv>0)));
% neighborAll.ArelFront = sqrt( (neighborAll.ArelFrontx)^2 + (neighborAll.ArelFronty)^2 + (neighborAll.ArelFrontz)^2);
% 
% neighborAll.ArelBackx = part0.ax - mean(tracer0.ax(idx01(rv<0)));
% neighborAll.ArelBacky = part0.ay - mean(tracer0.ay(idx01(rv<0)));
% neighborAll.ArelBackz = part0.az - mean(tracer0.az(idx01(rv<0)));
% neighborAll.ArelBack = sqrt( (neighborAll.ArelBackx)^2 + (neighborAll.ArelBacky)^2 + (neighborAll.ArelBackz)^2);


%% layer meshing of neighbors

R = Rmin:dR:Rmax; 

for k = 1:numel(R)-1
    neighborLayer(k).Rmin = Rmin;
    neighborLayer(k).Rmax = Rmax;
    neighborLayer(k).dR = dR;

    rlow = max(Rmin,R(k));
    rup  = min(Rmax,R(k+1));
    neighborLayer(k).rlow = rlow;
    neighborLayer(k).rup = rup;

    idx02 = find(d>rlow & d<rup);

    neighborLayer(k).idx = idx2(idx02);
    neighborLayer(k).d = d(idx02);
    
    clearvars rpt2 rv2 
    rpt2.x = tracer0.xf(idx02) - part0.xf;
    rpt2.y = tracer0.yf(idx02) - part0.yf;
    rpt2.z = tracer0.zf(idx02) - part0.zf;
    rv2 = part0.vx*rpt2.x + part0.vy*rpt2.y + part0.vz*rpt2.z;
    neighborLayer(k).idxfront = neighborLayer(k).idx(rv2>0);
    neighborLayer(k).idxback = neighborLayer(k).idx(rv2<0);

%     neighborLayer(k).VrelFrontx = part0.vx - mean(tracer0.vx((idx02(rv2>0))));
%     neighborLayer(k).VrelFronty = part0.vy - mean(tracer0.vy((idx02(rv2>0))));
%     neighborLayer(k).VrelFrontz = part0.vz - mean(tracer0.vz((idx02(rv2>0))));
%     neighborLayer(k).VrelFront = sqrt( (neighborLayer(k).VrelFrontx)^2 + (neighborLayer(k).VrelFronty)^2 + (neighborLayer(k).VrelFrontz)^2);
%    
%     neighborLayer(k).VrelBackx = part0.vx - mean(tracer0.vx((idx02(rv2<0))));
%     neighborLayer(k).VrelBacky = part0.vy - mean(tracer0.vy((idx02(rv2<0))));
%     neighborLayer(k).VrelBackz = part0.vz - mean(tracer0.vz((idx02(rv2<0))));
%     neighborLayer(k).VrelBack = sqrt( (neighborLayer(k).VrelBackx)^2 + (neighborLayer(k).VrelBacky)^2 + (neighborLayer(k).VrelBackz)^2);
%     
%     neighborLayer(k).ArelFrontx = part0.ax - mean(tracer0.ax((idx02(rv2>0))));
%     neighborLayer(k).ArelFronty = part0.ay - mean(tracer0.ay((idx02(rv2>0))));
%     neighborLayer(k).ArelFrontz = part0.az - mean(tracer0.az((idx02(rv2>0))));
%     neighborLayer(k).ArelFront = sqrt( (neighborLayer(k).ArelFrontx)^2 + (neighborLayer(k).ArelFronty)^2 + (neighborLayer(k).ArelFrontz)^2);
%     
%     neighborLayer(k).ArelBackx = part0.ax - mean(tracer0.ax((idx02(rv2<0))));
%     neighborLayer(k).ArelBacky = part0.ay - mean(tracer0.ay((idx02(rv2<0))));
%     neighborLayer(k).ArelBackz = part0.az - mean(tracer0.az((idx02(rv2<0))));
%     neighborLayer(k).ArelBack = sqrt( (neighborLayer(k).ArelBackx)^2 + (neighborLayer(k).ArelBacky)^2 + (neighborLayer(k).ArelBackz)^2);
end





