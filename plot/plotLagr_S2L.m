function plotLagr_S2L(S2L,Fs,mycolormap,ifsave,fout)

color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 2;

% figure;loglog(S2Lx.tau,S2Lx.mean./S2Lx.tau/Fs/2)
figure;
loglog(S2L(1).tau/Fs,S2L(1).mean,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=lw);hold on
loglog(S2L(2).tau/Fs,S2L(2).mean,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=lw);
loglog(S2L(3).tau/Fs,S2L(3).mean,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=lw);

lim1 = max(S2L(1).tau)/10;
lim2 = max(S2L(1).tau)/2;

xS2L = linspace(1,lim1,100)/Fs;
am1 = 20*min(ylim)/(min(xlim)^2);
loglog(xS2L,am1*xS2L.^2,'--',Color=color1,LineWidth=lw)

xS2L = linspace(lim1,lim2,100)/Fs;
am2 = am1*(lim1/Fs)^2/(lim1/Fs);
loglog(xS2L,am2*xS2L.^1,'--',Color=color1,LineWidth=lw)
% xS2L = linspace(100,300,100);
% loglog(xS2L,8e4*xS2L.^0,'--',Color=color1,LineWidth=lw)

set(gca,FontSize=15)
legend('$S_2^L(x)$','$S_2^L(y)$','$S_2^L(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^L$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^L$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau/s$','interpreter','latex',FontWeight='bold',FontSize=24)
lc1 = lim1/2/Fs;
lc2 = lim2/2/Fs;
text(lc1,2*am1*lc1^2,'$\tau^2$','interpreter','latex',FontWeight='bold',FontSize=18)
text(lc2,2*am2*lc2^1,'$\tau$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

if ifsave ==1
    figname = ['./' fout '/S2L'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end