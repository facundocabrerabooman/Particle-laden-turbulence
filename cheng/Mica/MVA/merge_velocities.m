function vel=merge_velocities(vel1,vel2);

%vel=merge_velocities(vel1,vel2);
if(~isfield(vel1.data,'acc'))
    for jj=1:numel(vel1.data)
        vel1.data(jj).acc=[];
        vel1.data(jj).velf=[];
    end
end
if(~isfield(vel2.data,'acc'))
    for jj=1:numel(vel2.data)
        vel2.data(jj).acc=[];
        vel2.data(jj).velf=[];
    end
end

vel.good=[vel1.good vel2.good+numel(vel1.data)];
vel.lmin=vel1.lmin;
vel.data(vel.good)=[vel1.data(vel1.good) vel2.data(vel2.good)];
vel.i0(vel.good)=[vel1.i0(vel1.good) vel2.i0(vel2.good)];
vel.length(vel.good)=[vel1.length(vel1.good) vel2.length(vel2.good)];
vel.file(vel.good)=[vel1.file(vel1.good) vel2.file(vel2.good)];
vel.status(vel.good)=[vel1.status(vel1.good) vel2.status(vel2.good)];