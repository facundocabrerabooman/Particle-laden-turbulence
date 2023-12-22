function plotEuler_PSDs(eulerStats,mycolormap,ifsave,fout)

color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 1;

figure
loglog(eulerStats.PSDk,eulerStats.PSD,'-',MarkerSize=8,Color=color3(1,:),LineWidth=lw);hold on

set(gca,FontSize=15)
% legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$PSD$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PSD$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$k$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

if ifsave ==1
    figname = ['./' fout '/PSD'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end