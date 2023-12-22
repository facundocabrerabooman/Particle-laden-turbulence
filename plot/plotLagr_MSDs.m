function plotLagr_MSDs(MSD,Fs,mycolormap,ifsave,fout)

color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 2;

figure;
loglog(MSD(1).tau/Fs,MSD(1).mean,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=lw);hold on
loglog(MSD(2).tau/Fs,MSD(2).mean,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=lw);
loglog(MSD(3).tau/Fs,MSD(3).mean,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=lw);


xMSD = linspace(1,max(MSD(1).tau),1000)/Fs;
loglog(xMSD,5*max(ylim)/(max(xlim)^2)*xMSD.^2,'--',Color=color1,LineWidth=lw)

set(gca,FontSize=15)
legend('MSDx','MSDy','MSDz',Location='best',FontSize=12)
title('$MSD$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$MSD(mm^2)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau(s)$','interpreter','latex',FontWeight='bold',FontSize=24)
text(max(xlim)/10,5*max(ylim)/(max(xlim)^2)*(max(xlim)/10)^2,'$\tau^2$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

if ifsave ==1
    figname = ['./' fout '/MSDs'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end