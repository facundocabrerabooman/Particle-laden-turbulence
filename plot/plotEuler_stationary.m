function plotEuler_stationary(eulerStats,Fs,mycolormap,ifsave,fout)

color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 1;

figure
t=tiledlayout(4,1,'TileSpacing','tight');
nexttile;
plot(eulerStats.Vmoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=lw);hold on
plot(eulerStats.VmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=lw);
plot(eulerStats.VmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=lw);
plot(eulerStats.VmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=lw);
set(gca,FontSize=15)
legend('$\langle \sqrt{x^2+y^2+z^2} \rangle$','$\langle x \rangle$','$\langle y \rangle$','$\langle z \rangle$','interpreter','latex',Location='best',FontSize=10);
title('$|V|, \sigma_V, |A|,\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$|V|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Vstd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=lw);hold on
plot(eulerStats.VstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=lw);
plot(eulerStats.VstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=lw);
plot(eulerStats.VstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=lw);
set(gca,FontSize=15)
ylabel('$\sigma_V$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Amoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=lw);hold on
plot(eulerStats.AmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=lw);
plot(eulerStats.AmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=lw);
plot(eulerStats.AmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=lw);
set(gca,FontSize=15)
ylabel('$|A|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Astd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=lw);hold on
plot(eulerStats.AstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=lw);
plot(eulerStats.AstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=lw);
plot(eulerStats.AstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=lw);
set(gca,FontSize=15)
ylabel('$\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/s$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticks(0:Fs/10:size(eulerStats.Astd,2))
xticklabels(num2cell([0:Fs/10:size(eulerStats.Astd,2)]/Fs))
grid on


linkaxes(t.Children,'x')

if ifsave == 1
    figname = ['./' fout '/VAt'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end