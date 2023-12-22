
function [sig_out, dsig_out, d2sig_out]=accelerations(sig_in,w,L)

% vel_out=accelerations(vel_in,w)
%
% calculates acceleration from velocities in structure vel_in
% the gaussian kernel is calculated from 
% kernel = velfiltcoef(filterwidth, filterlen)

N=floor(L/2);
Ns=numel(sig_in);

kerp=posfiltcoef(w,L);
kerv=velfiltcoef(w,L);
kera=accfiltcoef(w,L);

sig_out_tmp=conv(kerp,sig_in);
dsig_out_tmp=conv(kerv,sig_in);
d2sig_out_tmp=conv(kera,sig_in);

sig_out=sig_out_tmp(L:Ns);
dsig_out=dsig_out_tmp(L:Ns);
d2sig_out=d2sig_out_tmp(L:Ns);
end

