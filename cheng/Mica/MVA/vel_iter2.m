function [vel]=vel_iter2(vel,limit_Hess,amp_thresh,NbMic,K)

%[vel]=vel_iter2(vel,limit_Hess,amp_thresh,NbMic,K)

limit_Hess_iter=limit_Hess/50;
length_tresh=512;
vel=vel_trie(vel,limit_Hess,amp_thresh);

% dirname='.';
% sname=sprintf('%s/serie%i_vel_N7_K13_%0.5g.mat',dirname,serie,limit_Hess);
good=vel.good;
vel.good=[];
for j=1:length(good)
    [B,A]=butter(4,[min(vel.data(good(j)).freq(K:length(vel.data(good(j)).freq)-K)) max(vel.data(good(j)).freq(K:length(vel.data(good(j)).freq)-K))]*2);
    fdataf=filtfilt(B,A,vel.data(good(j)).seg);
    [freq,iHes]=MVA13mult1(conj(fdataf'),1,NbMic,K);

    vel.good=[vel.good good(j)];
    vel.data(good(j)).seg=vel.data(good(j)).seg;
    vel.data(good(j)).freq=freq;
    vel.data(good(j)).iHes=iHes;
    vel.i0(good(j))=vel.i0(good(j));
    vel.length(good(j))=length(freq);
    vel.file(good(j))=vel.file(good(j));
    vel.status(good(j))=1;
end
vel=vel_trie(vel,limit_Hess_iter,amp_thresh);
vel.param=[limit_Hess limit_Hess_iter amp_thresh NbMic K ];
ff=char(vel.file(1));
eval(['save vel_iter_' ff(1:6) '.mat vel;']);



