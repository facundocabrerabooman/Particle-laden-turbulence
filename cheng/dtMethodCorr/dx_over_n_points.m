function [dx] = dx_over_n_points(x,n)
%ACC_OVER_N_POINTS Summary of this function goes here
%   Detailed explanation goes here

%a(t)=x(t+n)+x(t-n)-2*x(t);
%a=x(2*n+1:end)+x(1:end-2*n)-2*x(n+1:end-n);

dx1=x(2*n+1:end)-x(1:end-2*n); % facteur 2 par rapoort au gradient
dx=dx1/2;
% dx2=x(2*(2*n)+1:end)-x(1:end-2*(2*n));
% dx=(8*dx1(n+1:end-n)-dx2)/12;


