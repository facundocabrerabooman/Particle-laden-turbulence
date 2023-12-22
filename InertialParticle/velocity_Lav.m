%%
tmin = 0.06;


%%
mycolormap = mycolor('#063970','#eeeee4','#e28743');%('#063970','#eeeee4','#e28743')
Nexp = size(part,1);
Nframemax = 1e6;
NTrackmax = 1e9;
mycolorind = 1:floor(size(mycolormap,1)/Nexp):size(mycolormap,1);
patchcolor = mycolormap(size(mycolormap,1)/2,:);
patchalpha = 0.4;

%%
figure
t=tiledlayout(3,1,'TileSpacing','none');
nexttile;
for nexp = 1:Nexp
    T = (part(nexp,:).T-(nexp-1)*Nframemax)./fps;
    semilogx(T,part(nexp,:).Vx,'-',Color=mycolormap(mycolorind(nexp),:),LineWidth=2);hold on
end
ylim([-0.35 0.35])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor)
alpha(patchalpha)
set(gca,FontSize=15)
title('$Velocity$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$V_x (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf 0.3])

%%
nexttile;
for nexp = 1:Nexp
    T = (part(nexp,:).T-(nexp-1)*Nframemax)./fps;
    semilogx(T,part(nexp,:).Vy,'-',Color=mycolormap(mycolorind(nexp),:),LineWidth=2);hold on
end
ylim([-0.7 0])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor)
alpha(patchalpha)
set(gca,FontSize=15)
ylabel('$V_y (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf 0.3])

%%
nexttile;
for nexp = 1:Nexp
    T = (part(nexp,:).T-(nexp-1)*Nframemax)./fps;
    semilogx(T,part(nexp,:).Vz,'-',Color=mycolormap(mycolorind(nexp),:),LineWidth=2);hold on
end
ylim([-0.35 0.35])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor)
alpha(patchalpha)
set(gca,FontSize=15)
ylabel('$V_z (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim([-inf 0.3])

%%
linkaxes(t.Children,'x')
colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.')
set(col,'position',[0.5 0.35 0.2 0.02])