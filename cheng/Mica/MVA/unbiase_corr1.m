function [R_out]=unbiase_corr1(R_in,sigs)

% R_out=unbiase_corr1(R_in,vel_in)
%
% - R_in is the correlation function non normalized to one at tau=zero
% (R_in(0)=sigma_tot)
% - vel_in is structure with velf (typically velm_xmm.data, or velm_xmm.data(velm_xdata..good))
%
% unbiases the correlation function R_in based on velocity statistics from
% vel_in following Nicolas Mordant's procedure :
%
% R_out(tau)=R_in(tau)/sigma_>(tau)
%
% where sigma_>(tau) is the variance of velocity signals longer than tau
%

I=(1:numel(R_in.mean));

R_out.mean=1-(sigs(1)-R_in.mean)./sigs(I);
R_out.std=R_in.std./sigs(I);
R_out.N=R_in.N;


	
	


