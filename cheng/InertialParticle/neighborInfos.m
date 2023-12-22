% function dist = relative_val(ppart,tpart,fields,Rmin,Rmax,Nframemax)
function neighbor = neighborInfos(part,tracer,idx1,idx3,Rmin,Rmax)

part0.x = part.Xf(idx1);
part0.y = part.Yf(idx1);
part0.z = part.Zf(idx1);
part0.vx = part.Vx(idx1);
part0.vy = part.Vy(idx1);
part0.vz = part.Vz(idx1);
part0.ax = part.Ax(idx1);
part0.ay = part.Ay(idx1);
part0.az = part.Az(idx1);

tracer0.x = tracer.Xf(idx3);
tracer0.y = tracer.Yf(idx3);
tracer0.z = tracer.Zf(idx3);
tracer0.vx = tracer.Vx(idx3);
tracer0.vy = tracer.Vy(idx3);
tracer0.vz = tracer.Vz(idx3);
tracer0.ax = tracer.Ax(idx3);
tracer0.ay = tracer.Ay(idx3);
tracer0.az = tracer.Az(idx3);

%% neighboring tracers

neighbor= struct('idx',[],'d',[],'Vrelx',[],'Vrely',[],'Vrelz',[],'Vrel',[],'Arelx',[],'Arely',[],'Arelz',[],'Arel',[]);

neighbor.Rmin = Rmin;
neighbor.Rmax = Rmax;

for i = 1:numel(tracer0.x)
    d(i,:) = sqrt((part0.x-tracer0.x(i))^2+(part0.y-tracer0.y(i))^2+(part0.z-tracer0.z(i))^2);
end
idx = find(d>Rmin & d<Rmax);
neighbor.idx = idx3(idx);
neighbor.d = d(idx);
neighbor.Vrelx = part0.vx - mean(tracer0.vx(idx));
neighbor.Vrely = part0.vy - mean(tracer0.vy(idx));
neighbor.Vrelz = part0.vz - mean(tracer0.vz(idx));
neighbor.Vrel = sqrt( (neighbor.Vrelx)^2 + (neighbor.Vrely)^2 + (neighbor.Vrelz)^2);
neighbor.Arelx = part0.ax - mean(tracer0.ax(idx));
neighbor.Arely = part0.ay - mean(tracer0.ay(idx));
neighbor.Arelz = part0.az - mean(tracer0.az(idx));
neighbor.Arel = sqrt( (neighbor.Arelx)^2 + (neighbor.Arely)^2 + (neighbor.Arelz)^2);
%% front and back 
rpt.x = tracer0.x(idx) - part0.x;
rpt.y = tracer0.y(idx) - part0.y;
rpt.z = tracer0.z(idx) - part0.z;

rv = part0.vx*rpt.x + part0.vy*rpt.y + part0.vz*rpt.z;
neighbor.idxfront = neighbor.idx(rv>0);
neighbor.idxback = neighbor.idx(rv<0);




