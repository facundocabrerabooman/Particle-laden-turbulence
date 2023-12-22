function [d2x] = d2x_over_n_points(x,n)
%ACC_OVER_N_POINTS Summary of this function goes here
%   Detailed explanation goes here


%a(t)=x(t+n)+x(t-n)-2*x(t);

d2x=x(2*n+1:end)+x(1:end-2*n)-2*x(n+1:end-n);

% 
% dx1=x(2*n+1:end)-x(1:end-2*n); % facteur 2 par rapoort au gradient
% dx2=x(2*(2*n)+1:end)-x(1:end-2*(2*n));
% dx=(8*df1(n+1:end-n)-df2)/12;

%d2fdt2=(f(2*n+1:end)+f(1:end-2*n)-2*f(1+n:end-n))/(n*dt)^2;

%d2x1=x(2*n+1:end)+x(1:end-2*n);%-2*f(1+n:end-n)
%d2x2=x(2*2*n+1:end)+x(1:end-2*2*n);%-2*f(1+n:end-n);
%d2x=(16*d2x1(n+1:end-n)-d2x2-30*x(1+2*n:end-2*n))/12;
