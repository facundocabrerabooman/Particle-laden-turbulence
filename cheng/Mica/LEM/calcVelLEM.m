function tracksout = calcVelLEM(tracksin,w,l,Fs)
%%
if exist('w','var')
    filterwidth = w;
    filterlength = l;
else
    filterwidth = 3;
    filterlength = 9;
end

if ~exist('Fs','var')
    Fs=1;
end

kerp = posfiltcoef(filterwidth,filterlength);
kerv = velfiltcoef(filterwidth,filterlength);
kera = accfiltcoef(filterwidth,filterlength);

%% Calc filtered position
tracksout = addStructFun(tracksin,'X','Xf',@(X)(conv(X,kerp,'valid')));
tracksout = addStructFun(tracksout,'Y','Yf',@(X)(conv(X,kerp,'valid')));
tracksout = addStructFun(tracksout,'Z','Zf',@(X)(conv(X,kerp,'valid')));
%% Calc filtered velocity
tracksout = addStructFun(tracksout,'X','Vx',@(X)(conv(X,kerv,'valid')*Fs));
tracksout = addStructFun(tracksout,'Y','Vy',@(X)(conv(X,kerv,'valid')*Fs));
tracksout = addStructFun(tracksout,'Z','Vz',@(X)(conv(X,kerv,'valid')*Fs));
%% Calc filtered acceleration
tracksout = addStructFun(tracksout,'X','Ax',@(X)(conv(X,kera,'valid')*Fs^2));
tracksout = addStructFun(tracksout,'Y','Ay',@(X)(conv(X,kera,'valid')*Fs^2));
tracksout = addStructFun(tracksout,'Z','Az',@(X)(conv(X,kera,'valid')*Fs^2));
%% Calc filtered times
if isfield(tracksout,'T')
    tracksout = addStructFun(tracksout,'T','Tf',@(X)(X(1+ceil((filterlength-1)/2):end-ceil((filterlength-1)/2))));
end
if isfield(tracksout,'Ntrack')
    tracksout = addStructFun(tracksout,'Ntrack','Ntrackf',@(X)(X(1+ceil((filterlength-1)/2):end-ceil((filterlength-1)/2))));
end
