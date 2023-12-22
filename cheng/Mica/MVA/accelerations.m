function vel_out=accelerations(vel_in,w)

% vel_out=accelerations(vel_in,w)
%
% calculates acceleration from velocities in structure vel_in
% the gaussian kernel is calculated from 
% kernel = velfiltcoef(filterwidth, filterlen)

kerv=posfiltcoef(w,5*w);
kera=velfiltcoef(w,5*w);


%l=length(kera);
vel_out=vel_in;
%vel_out.filtparam=[w,3*w];
for jj=1:numel(vel_in.data(vel_in.good))
   %disp(sprintf('%i',jj));
   %if isempty(vel_in.data(vel_in.good(jj)))
       acctmp=conv(kera,(vel_in.data(vel_in.good(jj)).freq));
       %acc=acctmp(l:length(acctmp)-l);
   
       veltmp=conv(kerv,(vel_in.data(vel_in.good(jj)).freq));
       %velf=veltmp(l:length(veltmp)-l);
		ll=length(vel_in.data(vel_in.good(jj)).freq);
       vel_out.data(vel_out.good(jj)).acc=acctmp(numel(kera):ll);
       vel_out.data(vel_out.good(jj)).velf=veltmp(numel(kerv):ll);
   
       %for kk=1:vel_in.length(jj)-numel(ker)
       %     acc(kk)=sum(ker'.*vel_in.data(jj).freq(kk:kk+numel(ker)-1));
       %     vel.data(jj).acc=acc;
       %end
   %end
end

