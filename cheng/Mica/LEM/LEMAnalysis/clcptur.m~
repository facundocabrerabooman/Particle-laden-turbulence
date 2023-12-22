dirParent='/Users/mlopez/Documents/exp/LEM/pospro/';
for k=2:2:12
dirname=[dirParent num2str(k)];
cd(dirname);
load datos.mat
%
nu=1e-6;
%%
sr23.xx=Sf.Sxx.sfabs(:,2)./((Sf.Sxx.r.*Sf.Sxx.scaler*1e-3).^(2/3))';
sr23.yy=Sf.Syy.sfabs(:,2)./((Sf.Syy.r.*Sf.Syy.scaler*1e-3).^(2/3))';
sr23.xy=Sf.Sxy.sfabs(:,2)./((Sf.Sxy.r.*Sf.Sxy.scaler*1e-3).^(2/3))';
sr23.yx=Sf.Syx.sfabs(:,2)./((Sf.Syx.r.*Sf.Syx.scaler*1e-3).^(2/3))';
figure; hold on;
plot (Sf.Sxx.r*Sf.Sxx.scaler, sr23.xx, 'k')
plot (Sf.Syy.r*Sf.Syy.scaler, sr23.yy, 'b')
plot (Sf.Sxy.r*Sf.Sxy.scaler, sr23.xy, 'r')
plot (Sf.Syx.r*Sf.Syx.scaler, sr23.yx, 'g')
%
mxsr23.xx=max(sr23.xx);
mxsr23.yy=max(sr23.yy);
mxsr23.xy=max(sr23.xy);
mxsr23.yx=max(sr23.yx);
%
eps.xx=(mxsr23.xx/2.1)^1.5;
eps.yy=(mxsr23.yy/2.1)^1.5;
eps.xy=(mxsr23.xy/2.1*3/4)^1.5;
eps.yx=(mxsr23.yx/2.1*3/4)^1.5;
%
eta.xx=(nu^3/eps.xx)^.25;
eta.yy=(nu^3/eps.yy)^.25;
eta.xy=(nu^3/eps.xy)^.25;
eta.yx=(nu^3/eps.yx)^.25;
%
tau_eta.xx=(nu/eps.xx)^.5;
tau_eta.yy=(nu/eps.yy)^.5;
tau_eta.xy=(nu/eps.xy)^.5;
tau_eta.yx=(nu/eps.yx)^.5;

%%
sr3.xx=-Sf.Sxx.sf(:,3)./((Sf.Sxx.r.*Sf.Sxx.scaler*1e-3))';
sr3.yy=-Sf.Syy.sf(:,3)./((Sf.Syy.r.*Sf.Syy.scaler*1e-3))';

figure;plot (Sf.Sxx.r*Sf.Sxx.scaler, sr3.xx,'b')
hold on;plot (Sf.Syy.r*Sf.Syy.scaler, sr3.yy,'r')
%%
mxsr3.xx=max(sr3.xx);
mxsr3.yy=max(sr3.yy);
%
eps3.xx=(mxsr3.xx*(5/4));
eps3.yy=(mxsr3.yy*(5/4));
%
eta3.xx=(nu^3/eps3.xx)^.25;
eta3.yy=(nu^3/eps3.yy)^.25;
%
tau_eta3.xx=(nu/eps3.xx)^.5;
tau_eta3.yy=(nu/eps3.yy)^.5;

%%
%figure;plot(Sf.Sxy.sfabs(:,2)./Sf.Sxx.sfabs(:,2));
%
save('result.mat','eps','eps3','eta','eta3','tau_eta','tau_eta3');
end
