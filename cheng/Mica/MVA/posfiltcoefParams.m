function [B,A] = posfiltcoefParams(filterwidth, filterlen)

N = filterlen;
L = floor(N/2);
w = filterwidth;
t = [-L:L]';
B=exp(-L^2/w^2)/(2*L*exp(-L^2/w^2)-w*sqrt(pi)*erf(L/w));
A=(1-2*B*L)/(w*sqrt(pi)*erf(L/w));