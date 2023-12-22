function [CD]=CD_Brown(Re)
%
%tau_p=Taup_Brown(Re,Gamma,d,nu)
% 
CD=24./Re.*(1+0.150*Re.^0.681)+0.407./(1+8710./Re);
