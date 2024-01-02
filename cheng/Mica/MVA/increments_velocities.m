function [incr]=increments_velocities(vel_in,dt,n,w)

%[incr]=increments_velocities(vel_in,dt,n,w)

%c=340;
%fact=c/2/nu0/sin(theta/2)*32768;

incr.dt=dt;

for kk=1:numel(dt)
    kk
    ig=0;
    lc=0;
    for jj=1:numel(vel_in.good)
        if(length(vel_in.data(vel_in.good(jj)).velf)-dt(kk)>1)
            ig=ig+1;
            ind=1:length(vel_in.data(vel_in.good(jj)).velf)-dt(kk);
            data(ig).x=vel_in.data(vel_in.good(jj)).velf(ind+dt(kk))-vel_in.data(vel_in.good(jj)).velf(ind);
            data(ig).l=numel(data(ig).x);
            lc=lc+data(ig).l;
            data(ig).lc=lc;
            %inc=[inc vel_in.data(vel_in.good(jj)).velf(ind+dt(kk))-vel_in.data(vel_in.good(jj)).velf(ind)];
        end
    end
    
    %pdf=mkpdf4(data,n,w);
    pdf=mkpdf4(data,n);
    %[incr(kk).xpdf,incr(kk).pdf]=mkpdf3(inc,65,-5,5);
    %stat=stat_velocities2(vel_in,65);
    incr.pdf(kk).xpdf=pdf.xpdf;
    incr.pdf(kk).pdf=pdf.pdf;
    incr.pdf(kk).xpdfn=pdf.xpdfn;
    incr.pdf(kk).pdfn=pdf.pdfn;
    incr.N(kk)=pdf.N;
    incr.mean(kk)=pdf.mean;
    incr.meanb(kk)=pdf.meanb;
    incr.std(kk)=pdf.std;
    incr.stdb(kk)=pdf.stdb;
    incr.flatness(kk)=pdf.flatness;
    incr.flatnessb(kk)=pdf.flatnessb;
    incr.skewness(kk)=pdf.skewness;
    incr.skewnessb(kk)=pdf.skewnessb;
    clear data;
    %eval(['vel_out.inc_' num2str(dt(kk)) '.x=x;']);
    %eval(['vel_out.inc_' num2str(dt(kk)) '.pdf=pdf;']);
end

