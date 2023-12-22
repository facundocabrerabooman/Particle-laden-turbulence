function plotLagr_std(w,s,wopt,colormap,ifsave,fout)

color1 = '#a155b9';
color3 = [colormap(1,:);colormap((size(colormap,1)+1)/2,:);colormap(end,:)];
lw= 2;

figure;
yyaxis left
loglog(w,s(1).vx,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=lw);hold on;
loglog(w,s(2).vx,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=lw);
loglog(w,s(3).vx,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=lw);
hold off

yyaxis right
loglog(w,s(1).ax,'o-',MarkerSize=8,Color=color3(1,:),LineWidth=lw);hold on;
loglog(w,s(2).ax,'o-',MarkerSize=8,Color=color3(2,:),LineWidth=lw);
loglog(w,s(3).ax,'o-',MarkerSize=8,Color=color3(3,:),LineWidth=lw);

plot([wopt wopt],ylim,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
yyaxis left
legend('$V_x$','$V_y$','$V_z$','$A_x$','$A_y$','$A_z$','interpreter','latex',Location='southwest',FontSize=12);
title('$std.(w)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$fliter\ width\ w$','interpreter','latex',FontWeight='bold',FontSize=18)

yyaxis left
ylabel('$\sigma_{v}$','interpreter','latex',FontWeight='bold',FontSize=24)
yyaxis right
ylabel('$\sigma_{a}$','interpreter','latex',FontWeight='bold',FontSize=24)

grid on
axis padded

if ifsave==1
    figname = ['./' fout '/Std'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end