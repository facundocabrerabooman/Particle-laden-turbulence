function stat=stat_particles(vel)
stat.mean=zeros(1,numel(vel.good));
stat.std=zeros(1,numel(vel.good));
stat.length=zeros(1,numel(vel.good));

for kk=1:numel(vel.good)
    stat.mean(kk)=mean(vel.data(vel.good(kk)).velf);
    stat.std(kk)=std(vel.data(vel.good(kk)).velf);
    stat.length(kk)=length(vel.data(vel.good(kk)).velf);
end