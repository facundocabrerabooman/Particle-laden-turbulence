function [S2batch,eta,L,Lint,Ceps] = batchelorS2(r,nu,sigmaU,epsilon,zeta)

%% S2batch = batchelorS2(r,a,L,espilon,nu,zeta)
%
% Batchelor parametrization for Eulerian 2nd ordre structure function, from
% https://www.researchgate.net/deref/http%3A%2F%2Farxiv.org%2Fabs%2Fchao-dyn%2F9704015v2?_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6Il9kaXJlY3QiLCJwYWdlIjoicHVibGljYXRpb24iLCJwb3NpdGlvbiI6InBhZ2VDb250ZW50In19
%
% use zeta = 2/3 for K41
% a and L are fit parameters 

%L = (sigmaU^2*30*nu/epsilon*a^(zeta-2))^(1/zeta);
%S2batch = epsilon/15/nu * r.^2./(1+r.^2/a^2).^(1-zeta/2).*1./(1+r.^2/L^2).^(zeta/2);

eta = (nu^3/epsilon).^0.25
C2 = 2.1
L = (2*sigmaU^2/C2)^(1/zeta)/epsilon;
a = (15/C2)^(1/(2-zeta))*eta;
% Lint = Ceps*L;


S2batch = C2*epsilon^zeta/a^(2-zeta) * r.^2./(1+r.^2/a^2).^(1-zeta/2).*1./(1+r.^2/L^2).^(zeta/2);

%%
rfit = linspace(0,1,1024);
Rint = cumsum(1-C2*epsilon^zeta/a^(2-zeta) * rfit.^2./(1+rfit.^2/a^2).^(1-zeta/2).*1./(1+rfit.^2/L^2).^(zeta/2)/2/sigmaU^2);
Lint = Rint(end)*rfit(2);
Ceps = Lint/L;




