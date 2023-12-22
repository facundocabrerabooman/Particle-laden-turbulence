function plotEuler_Sau(eulerStats,mycolormap,ifsave,fout)

color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 2;


figure
semilogx(eulerStats.r,eulerStats.Sau,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=lw);hold on
semilogx(eulerStats.r,eulerStats.Saulong,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=lw);hold on

set(gca,FontSize=15)
legend('$S_{au}$','$S_{au}^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

if ifsave ==1
    figname = ['./' fout '/Sau'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end