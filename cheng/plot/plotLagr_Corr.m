function plotLagr_Corr(Ruu,Raa,Fs,n,if2layers,mycolormap,ifsave,fout)

color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 1;


%% Correlation fit - Thomas's Script
Ruufit(1) = correlationFit(Ruu(1),Fs,1,n(1),'V',if2layers(1),1);
Ruufit(2) = correlationFit(Ruu(2),Fs,1,n(1),'V',if2layers(1),0);
Ruufit(3) = correlationFit(Ruu(3),Fs,1,n(1),'V',if2layers(1),1);

Raafit(1) = correlationFit(Raa(1),Fs,1,n(2),'A',if2layers(2),1);
Raafit(2) = correlationFit(Raa(2),Fs,1,n(2),'A',if2layers(2),1);
Raafit(3) = correlationFit(Raa(3),Fs,1,n(2),'A',if2layers(2),0);

%% plots
% figure;
% main plot: zoom in
f1 = figure;
tiledlayout(2,1)
nexttile
plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',MarkerSize=6,Color=color3(1,:),LineWidth=lw);hold on
plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',MarkerSize=6,Color=color3(2,:),LineWidth=lw);
plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',MarkerSize=6,Color=color3(3,:),LineWidth=lw);

plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:),LineWidth=2*lw);hold on
plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:),LineWidth=2*lw)
plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:),LineWidth=2*lw)

set(gca,FontSize=15)
legend('$R_{uu}(x)$','$R_{uu}(y)$','$R_{uu}(z)$','interpreter','latex',Location='best',FontSize=12)
title('$R_{uu}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$R_{uu}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau$/s','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis tight
% xlim([0 0.045])
% ylim([-1 inf])

nexttile
plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'o',MarkerSize=6,Color=color3(1,:),LineWidth=lw);hold on
plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'o',MarkerSize=6,Color=color3(2,:),LineWidth=lw);
plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'o',MarkerSize=6,Color=color3(3,:),LineWidth=lw);

plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:),LineWidth=2*lw)
plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:),LineWidth=2*lw)
plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:),LineWidth=2*lw)


set(gca,FontSize=15)
legend('$R_{aa}(x)$','$R_{aa}(y)$','$R_{aa}(z)$','interpreter','latex',Location='best',FontSize=12)
title('$R_{aa}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$R_{aa}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau$/s','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis tight
% xlim([0 0.045])
% ylim([-1 inf])


% add inset: zoom out 
% axes('Position',[0.4 0.3 0.3 0.2]);
% plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',MarkerSize=3,Color=color3(1,:),LineWidth=lw);hold on
% plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',MarkerSize=3,Color=color3(2,:),LineWidth=lw);
% plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',MarkerSize=3,Color=color3(3,:),LineWidth=lw);
% 
% plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'^',MarkerSize=3,Color=color3(1,:),LineWidth=lw);
% plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'^',MarkerSize=3,Color=color3(2,:),LineWidth=lw);
% plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'^',MarkerSize=3,Color=color3(3,:),LineWidth=lw);
% 
% plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:),LineWidth=lw);hold on
% plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:),LineWidth=lw)
% plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:),LineWidth=lw)
% 
% plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:),LineWidth=lw)
% plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:),LineWidth=lw)
% plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:),LineWidth=lw)
% set(gca,FontSize=12)
% grid on
% axis tight
% 
% ylim([-1 inf])

if ifsave ==1
    figname = ['./' fout '/Corr'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end