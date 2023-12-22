function [Ruu,Raa] = LagrangianStatCorr(track)

Ruu(1) = xcorr_struct(track,'Vx',1);
Ruu(2) = xcorr_struct(track,'Vy',1);
Ruu(3) = xcorr_struct(track,'Vz',1);

Raa(1) = xcorr_struct(track,'Ax',1);
Raa(2) = xcorr_struct(track,'Ay',1);
Raa(3) = xcorr_struct(track,'Az',1);