function [freq,iHes]=vel_iter_test(seg,freq,NbMic,K)

%[vel]=vel_iter2(vel,limit_Hess,amp_thresh,NbMic,K)

% dirname='.';
% sname=sprintf('%s/serie%i_vel_N7_K13_%0.5g.mat',dirname,serie,limit_Hess);

    [B,A]=butter(4,[1.*min(freq(K:length(freq)-K)) 1.*max(freq(K:length(freq)-K))]*2);
    fdataf=filtfilt(B,A,seg);
    [freq,iHes]=MVA13mult1(fdataf,1,NbMic,K);




