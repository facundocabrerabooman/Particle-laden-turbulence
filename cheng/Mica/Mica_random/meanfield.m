[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(track,'X',[2 4 6 8 ],[20 21 22],1,1,'x','cart');
[mYdt, mBxdt, bins, N] = track2meanDxDt3DProfile(track,'Y',[2 4 6 8 ],[20 21 22],1,1,'y','cart');
[mZdt, mBxdt, bins, N] = track2meanDxDt3DProfile(track,'Z',[2 4 6 8 ],[20 21 22],1,1,'z','cart');

I = find(N<20);

% mXdt(I)=NaN;
% mYdt(I)=NaN;
% mZdt(I)=NaN;

mXdts=smoothn(mXdt,0.1);
mYdts=smoothn(mYdt,0.1);
mZdts=smoothn(mZdt,0.1);

mXdts(I)=NaN;
mYdts(I)=NaN;
mZdts(I)=NaN;

%%
[X,Y,Z]=ndgrid(bins{1},bins{2},bins{3});

Fx=griddedInterpolant(X,Y,Z,mXdts,'linear','none');
Fy=griddedInterpolant(X,Y,Z,mYdts,'linear','none');
Fz=griddedInterpolant(X,Y,Z,mZdts,'linear','none');

%%
xx=linspace(-40,40,30);
yy=linspace(-40,40,30);
zz=linspace(-40,40,30);

[XX,YY,ZZ]=ndgrid(xx,yy,zz);

mVx = Fx(XX,YY,ZZ);
mVy = Fy(XX,YY,ZZ);
mVz = Fz(XX,YY,ZZ);

%%

[XX,YY,ZZ]=meshgrid(yy,xx,zz);
figure;slice(XX,YY,ZZ,mVx,0,0,0);