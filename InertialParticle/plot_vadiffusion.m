function plot_vadiffusion(part,axisrange,fps,ifsave,fout)

% axisrange = [0.05 0.3];

%%
Nframemax = 1e6;
NTrackmax = 1e9;
for i = 1:numel(part)
    nexp(i,:) = fix(part(i).T/Nframemax)+1;
    t(i,:) = mod(part(i).T,Nframemax)/fps;
end
uni_exp = unique(nexp);
Nexp = numel(uni_exp);

%%
tmin = min(t);
tmax = max(t);
tinterp = linspace(tmin,tmax,1000);

%%
mycolormap = mycolor('#063970','#eeeee4','#e28743');%('#063970','#eeeee4','#e28743')
mycolorind = floor(linspace(1,size(mycolormap,1),Nexp));
% patchcolor1 = mycolormap(size(mycolormap,1)/2,:);
patchcolor2 = mycolor('#a155b9');
% patchalpha1 = 0.4;
patchalpha2 = 0.3;

%% Velocity
Vinterp = struct('x',[],'y',[],'z',[]);

figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% vx
nexttile;
for i = 1:Nexp
    ind = find(nexp == uni_exp(i));
    Vinterp.x = [Vinterp.x; interp1(t(ind),vertcat(part(ind).Vx),tinterp)];
    plot(t(ind),vertcat(part(ind).Vx),'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.x,"omitnan") + std(Vinterp.x,"omitnan");
lowbound = mean(Vinterp.x,"omitnan") - std(Vinterp.x,"omitnan");
hm = plot(tinterp,mean(Vinterp.x,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-0.15 0.15])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
set(gca,FontSize=15)
legend([hm,hs],'$\langle V_i \rangle$','$\langle V_i \rangle \pm \sigma (V_i)$','interpreter','latex',Location='eastoutside',FontSize=12);
title('$Velocity$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$V_x (mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
grid on
xlim(axisrange)

% vy
nexttile;
for i = 1:Nexp
    ind = find(nexp == uni_exp(i));
    Vinterp.y = [Vinterp.y; interp1(t(ind),vertcat(part(ind).Vy),tinterp)];
    plot(t(ind),vertcat(part(ind).Vy),'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.y,"omitnan") + std(Vinterp.y,"omitnan");
lowbound = mean(Vinterp.y,"omitnan") - std(Vinterp.y,"omitnan");
hm = plot(tinterp,mean(Vinterp.y,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-0.7 0])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
set(gca,FontSize=15)
ylabel('$V_y (mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
grid on
xlim(axisrange)

% vz
nexttile;
for i = 1:Nexp
    ind = find(nexp == uni_exp(i));
    Vinterp.z = [Vinterp.z; interp1(t(ind),vertcat(part(ind).Vz),tinterp)];
    plot(t(ind),vertcat(part(ind).Vz),'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.z,"omitnan") + std(Vinterp.z,"omitnan");
lowbound = mean(Vinterp.z,"omitnan") - std(Vinterp.z,"omitnan");
hm = plot(tinterp,mean(Vinterp.z,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-0.15 0.15])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
set(gca,FontSize=15)
ylabel('$V_z (mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='eastoutside');
ylabel(col,'Num. Exp.')
% set(col,'position',[0.42 0.34 0.2 0.02])

%
% linkaxes(tf.Children,'x')

if ifsave ==1
    figname1 = [fout '/V'];
    savefig(figname1)
    saveas(gcf,figname1,'png')
    saveas(gcf,figname1,'pdf')
end
%% Acceleration
Ainterp = struct('x',[],'y',[],'z',[]);

figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% ax
nexttile;
for i = 1:Nexp
    ind = find(nexp == uni_exp(i));
    Ainterp.x = [Ainterp.x; interp1(t(ind),vertcat(part(ind).Ax),tinterp)];
    plot(t(ind),vertcat(part(ind).Ax),'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.x,"omitnan") + std(Ainterp.x,"omitnan");
lowbound = mean(Ainterp.x,"omitnan") - std(Ainterp.x,"omitnan");
hm = plot(tinterp,mean(Ainterp.x,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-250 250])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
set(gca,FontSize=15)
legend([hm,hs],'$\langle A_i \rangle$','$\langle A_i \rangle \pm \sigma (A_i)$','interpreter','latex',Location='eastoutside',FontSize=12);
title('$Acce.$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$A_x (mm/s^{2})$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
grid on
xlim(axisrange)

% ay
nexttile;
for i = 1:Nexp
    ind = find(nexp == uni_exp(i));
    Ainterp.y = [Ainterp.y; interp1(t(ind),vertcat(part(ind).Ay),tinterp)];
    plot(t(ind),vertcat(part(ind).Ay),'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.y,"omitnan") + std(Ainterp.y,"omitnan");
lowbound = mean(Ainterp.y,"omitnan") - std(Ainterp.y,"omitnan");
hm = plot(tinterp,mean(Ainterp.y,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-250 250])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
set(gca,FontSize=15)
ylabel('$A_y (mm/s^{2})$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
grid on
xlim(axisrange)

% az
nexttile;
for i = 1:Nexp
    ind = find(nexp == uni_exp(i));
    Ainterp.z = [Ainterp.z; interp1(t(ind),vertcat(part(ind).Az),tinterp)];
    plot(t(ind),vertcat(part(ind).Az),'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.z,"omitnan") + std(Ainterp.z,"omitnan");
lowbound = mean(Ainterp.z,"omitnan") - std(Ainterp.z,"omitnan");
hm = plot(tinterp,mean(Ainterp.z,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-250 250])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
set(gca,FontSize=15)
ylabel('$A_z (mm/s^{2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='eastoutside');
ylabel(col,'Num. Exp.');
% set(col,'position',[0.42 0.34 0.2 0.02])

%
% linkaxes(tf.Children,'x')

if ifsave ==1
    figname2 = [fout '/A'];
    savefig(figname2)
    saveas(gcf,figname2,'png')
    saveas(gcf,figname2,'pdf')
end

%% estimate taup
vytau = mean(Vinterp.y,"omitnan");
figure;
subplot(2,1,1)
semilogx(tinterp,vytau,'k.')
set(gca,FontSize=15);grid
ylabel('$\langle V_y(\tau) \rangle $','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau(s)$','interpreter','latex',FontWeight='bold',FontSize=18)

v0 = vytau(1);
vs = mean(vytau(tinterp<0.22 & tinterp>0.16));
vn = (vytau-v0)/(vs-v0);
subplot(2,1,2);
semilogx(tinterp,vn,'k.')
set(gca,FontSize=15);grid
ylabel('$(\langle V_y \rangle (\tau)-V_0)/(V_s-V_0)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
if ifsave ==1
    figname3 = [fout '/Vs'];
    savefig(figname3)
    saveas(gcf,figname3,'png')
    saveas(gcf,figname3,'pdf')
end
