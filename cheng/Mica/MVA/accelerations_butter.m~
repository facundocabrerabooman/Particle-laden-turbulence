function vel_out=accelerations_butter(vel_in,Fc,F_Ech)

% calculates acceleration from velocities in structure vel_in
% the gaussian kernel is calculated from 
% kernel = velfiltcoef(filterwidth, filterlen)

%kerv=posfiltcoef(w,3*w);
kera=velfiltcoef(w,3*w);

[B,A]=butter(4,Fc*2/F_Ech);

%l=length(kera);
vel_out=vel_in;
for jj=1:numel(vel_in.data)
   disp(sprintf('%i',jj));

   veltmp=filfilt(vel_in.data(jj).freq,B,A);
   %veltmp=conv(kerv,vel_in.data(jj).freq);
   %velf=veltmp(l:length(veltmp)-l);

   acctmp=conv(kera,vel_in.data(jj).freq);
   %acc=acctmp(l:length(acctmp)-l);
   
   
   vel_out.data(jj).acc=acctmp(numel(kera):vel_in.length(jj));
   vel_out.data(jj).velf=veltmp(numel(kerv):vel_in.length(jj));
   
   %for kk=1:vel_in.length(jj)-numel(ker)
   %     acc(kk)=sum(ker'.*vel_in.data(jj).freq(kk:kk+numel(ker)-1));
   %     vel.data(jj).acc=acc;
   %end
end

