function plotEuler_Spn(eulerStats,mycolormap,ifsave,fout)

color1 = '#a155b9';
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];
lw= 1;

figure;
loglog(eulerStats.r,eulerStats.Splong{1,1},'d-',MarkerSize=8,Color=color5(1,:),LineWidth=lw);hold on
loglog(eulerStats.r,eulerStats.Splong{1,2},'d-',MarkerSize=8,Color=color5(2,:),LineWidth=lw);
loglog(eulerStats.r,eulerStats.Splong{1,3},'d-',MarkerSize=8,Color=color5(3,:),LineWidth=lw);
loglog(eulerStats.r,eulerStats.Splong{1,4},'d-',MarkerSize=8,Color=color5(4,:),LineWidth=lw);
loglog(eulerStats.r,eulerStats.Splong{1,5},'d-',MarkerSize=8,Color=color5(5,:),LineWidth=lw);

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1,LineWidth=lw)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1,LineWidth=lw)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=lw)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=lw)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1,LineWidth=lw)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

if ifsave ==1
    figname = ['./' fout '/Splong'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end