function vel_out=split(vel)

for jj=1:numel(vel.good)
    fname=char(vel.file(vel.good(jj)));
    jjname(jj)=str2num(fname(8:numel(fname)));
end

[d,k]=find(diff(jjname)<0);

k=[0  k  numel(jjname)];

for jj=1:numel(k)-1 
    vel_out(jj).lmin=vel.lmin;
    vel_out(jj).good=1:k(jj+1)-k(jj);
    vel_out(jj).data=vel.data(vel.good(k(jj)+1:k(jj+1)));
    vel_out(jj).i0=vel.i0(vel.good(k(jj)+1:k(jj+1)));
    vel_out(jj).length=vel.length(vel.good(k(jj)+1:k(jj+1)));
    vel_out(jj).file=vel.file(vel.good(k(jj)+1:k(jj+1)));
    vel_out(jj).status=vel.status(vel.good(k(jj)+1:k(jj+1)));
end