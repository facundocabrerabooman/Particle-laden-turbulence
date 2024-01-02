function stat=stat_velocities(vel)

data=[]

for j=1:length(vel)
    data=[data vel(j).freq];
end

stat.mean=mean(data);
stat.std=std(data);
stat.hist=hist(data,64);