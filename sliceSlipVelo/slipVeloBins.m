function [mygrids,slipVelo] = slipVeloBins(track,Fs,dt,nbins,n,power)

%% mean fields
[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(track,'Xf',dt,nbins,n,power,'x','cart');
[mYdt, mBydt, ~, ~] = track2meanDxDt3DProfile(track,'Yf',dt,nbins,n,power,'y','cart');
[mZdt, mBzdt, ~, ~] = track2meanDxDt3DProfile(track,'Zf',dt,nbins,n,power,'z','cart');

% as dt is intergra array (like [2 4 6 8]), we have to multiply it by Fs^n
mXdt=mXdt*Fs^n;
mYdt=mYdt*Fs^n;
mZdt=mZdt*Fs^n;


% threshold = 10;
% I = find(N<threshold);

% mXdts=smoothn(mXdt,0.1);
% mYdts=smoothn(mYdt,0.1);
% mZdts=smoothn(mZdt,0.1);
% 
% mXdts(I)=NaN;
% mYdts(I)=NaN;
% mZdts(I)=NaN;

%% prepare for slice visualization
[X,Y,Z]=ndgrid(bins{1},bins{2},bins{3});
[mygrids.XX,mygrids.YY,mygrids.ZZ]=meshgrid(Y,X,Z); 

slipVelo.XX = mXdt;
slipVelo.YY = mYdt;
slipVelo.ZZ = mZdt;
