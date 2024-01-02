function vel_out=increments(vel_in,dt,nu0,theta)

c=340;
fact=c/2/nu0/sin(theta/2)*32768;
vel_out=vel_in;

for kk=1:numel(dt)
    inc=[];
    for jj=1:numel(vel_in.good)

        ind=1:vel_in.length(vel_in.good(jj))-dt(kk);
        inc=[inc vel_in.data(vel_in.good(jj)).freq(ind+dt)-vel_in.data(vel_in.good(jj)).freq(ind)];
    end
    [x,pdf]=mkpdf3(inc*fact,32,-5,5);
    eval(['vel_out.inc_' num2str(dt(kk)) '.x=x;']);
    eval(['vel_out.inc_' num2str(dt(kk)) '.pdf=pdf;']);
end

