function [mygrids,mymeanFields,trackout] = meanFields(track,Fs,dt,nbins,threshold,n,power,gridRange,ifsave)

%% mean fields
[mXdt, mBxdt, bins, N] = track2meanDxDt3DProfile(track,'Xf',dt,nbins,n,power,'x','cart');
[mYdt, mBydt, ~, ~] = track2meanDxDt3DProfile(track,'Yf',dt,nbins,n,power,'y','cart');
[mZdt, mBzdt, ~, ~] = track2meanDxDt3DProfile(track,'Zf',dt,nbins,n,power,'z','cart');

% as dt is intergra array (like [2 4 6 8]), we have to multiply it by Fs^n
mXdt=mXdt*Fs^n;
mYdt=mYdt*Fs^n;
mZdt=mZdt*Fs^n;


% threshold = 10;
I = find(N<threshold);

mXdts=smoothn(mXdt,0.1);
mYdts=smoothn(mYdt,0.1);
mZdts=smoothn(mZdt,0.1);

mXdts(I)=NaN;
mYdts(I)=NaN;
mZdts(I)=NaN;

%% interpolation
[X,Y,Z]=ndgrid(bins{1},bins{2},bins{3});

Fx=griddedInterpolant(X,Y,Z,mXdts,'linear','none');
Fy=griddedInterpolant(X,Y,Z,mYdts,'linear','none');
Fz=griddedInterpolant(X,Y,Z,mZdts,'linear','none');

xx=linspace(gridRange.x(1),gridRange.x(2),100);
yy=linspace(gridRange.y(1),gridRange.y(2),200);
zz=linspace(gridRange.z(1),gridRange.z(2),100);

[XX,YY,ZZ]=ndgrid(xx,yy,zz);

mymeanFields.x = Fx(XX,YY,ZZ);
mymeanFields.y = Fy(XX,YY,ZZ);
mymeanFields.z = Fz(XX,YY,ZZ);

%% prepare for slice visualization
[mygrids.XX,mygrids.YY,mygrids.ZZ]=meshgrid(yy,xx,zz); 

%%
if ifsave ==1
    disp('saving')
    if n==1 
        if power ==1
            filename = './Meanfields_Vel.mat'
        elseif power ==2
            filename = './Meanfields_Vel_rms.mat'
        end
    elseif n==2
        if power ==1
            filename = './Meanfields_Acce.mat'
        elseif power ==2
            filename = './Meanfields_Acce_rms.mat'
        end
    end
    
    save(filename,'mXdt','mYdt','mZdt','mBxdt','mBydt','mBzdt','bins','N','threshold','mXdts','mYdts','mZdts','mymeanFields','mygrids')
end
%% substact mean fields
trackout = [];
if n==1 && power ==1
    trackout = substractMeanField(track,Fx,'Vx',{'mVx','sVx'});
    trackout = substractMeanField(trackout,Fy,'Vy',{'mVy','sVy'});
    trackout = substractMeanField(trackout,Fz,'Vz',{'mVz','sVz'});
elseif n==2 && power ==1
    trackout = substractMeanField(track,Fx,'Ax',{'mAx','sAx'});
    trackout = substractMeanField(trackout,Fy,'Ay',{'mAy','sAy'});
    trackout = substractMeanField(trackout,Fz,'Az',{'mAz','sAz'});
else
    trackout = track;
end