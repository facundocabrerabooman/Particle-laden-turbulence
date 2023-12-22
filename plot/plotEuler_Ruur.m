function plotEuler_Ruur(eulerStats,mycolormap,ifsave,fout)

color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 1;

figure;
plot(eulerStats.Ruur,eulerStats.Ruu,'-',MarkerSize=8,Color=color3(1,:),LineWidth=lw);hold on
set(gca,FontSize=15)
% legend('$Ruu(r)$','interpreter','latex',Location='best',FontSize=12)
title('$Ruu(r)$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$Ruu(r)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$Ruu_r$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

if ifsave ==1
    figname = ['./' fout '/Ruu'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end