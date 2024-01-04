function tracksout = calcVelLEM(tracksin,w,l,Fs)
%%
if exist('w','var')
    filterwidth = w;
    filterlength = l;

else
%     filterwidth = 3;
%     filterlength = 9;
disp('Set filter width')
end

if ~exist('Fs','var')
    Fs=1;
end

kerp = posfiltcoef(filterwidth,filterlength);
kerv = velfiltcoef(filterwidth,filterlength);
kera = accfiltcoef(filterwidth,filterlength);
%% Calc filtered position
tracksout = addStructFun(tracksin,'x','Xf',@(X)(conv(X,kerp,'valid')));
tracksout = addStructFun(tracksout,'y','Yf',@(X)(conv(X,kerp,'valid')));
tracksout = addStructFun(tracksout,'z','Zf',@(X)(conv(X,kerp,'valid')));
%% Calc filtered velocity
tracksout = addStructFun(tracksout,'x','Vx',@(X)(conv(X,kerv,'valid')*Fs));
tracksout = addStructFun(tracksout,'y','Vy',@(X)(conv(X,kerv,'valid')*Fs));
tracksout = addStructFun(tracksout,'z','Vz',@(X)(conv(X,kerv,'valid')*Fs));
%% Calc filtered acceleration
tracksout = addStructFun(tracksout,'x','Ax',@(X)(conv(X,kera,'valid')*Fs^2));
tracksout = addStructFun(tracksout,'y','Ay',@(X)(conv(X,kera,'valid')*Fs^2));
tracksout = addStructFun(tracksout,'z','Az',@(X)(conv(X,kera,'valid')*Fs^2));
%% Calc filtered times
if isfield(tracksout,'t')
    tracksout = addStructFun(tracksout,'t','Tf',@(X)(X(1+ceil((filterlength-1)/2):end-ceil((filterlength-1)/2))));
    tracksout = addStructFun(tracksout,'t','Tf_acc',@(X)(X(1+ceil((filterlength-1)/2):end-ceil((filterlength-1)/2))));
end
if isfield(tracksout,'Ntrack')
    tracksout = addStructFun(tracksout,'Ntrack','Ntrackf',@(X)(X(1+ceil((filterlength-1)/2):end-ceil((filterlength-1)/2))));
    tracksout = addStructFun(tracksout,'Ntrack','Ntrack_acc',@(X)(X(1+ceil((filterlength-1)/2):end-ceil((filterlength-1)/2))));
end
end
