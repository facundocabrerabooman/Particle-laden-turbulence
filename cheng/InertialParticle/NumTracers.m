clear
clc
close all

%%
fps = 4230;

ppart = load('.\particle\STB.mat');
tpart = load('.\tracers\EulerianPart.mat');

%%
Nframemax = 1e6;
pNexp = numel(unique(fix(ppart.part.T/Nframemax)+1));

%%
mycolormap = mycolor('#063970','#eeeee4','#e28743');%('#063970','#eeeee4','#e28743')
mycolorind = 1:floor(size(mycolormap,1)/pNexp):size(mycolormap,1);
patchcolor1 = mycolormap(size(mycolormap,1)/2,:);
patchcolor2 = mycolor('#a155b9');
patchalpha1 = 0.4;
patchalpha2 = 0.3;

%%
dist = relative_val(ppart,tpart,'all',0,5,Nframemax);


mean_num  = mean([dist.num]);
mean_d  = mean(arrayfun(@(X)(mean(vertcat(X.d{1,:}))),dist));
disp(['Mean num   =  ' num2str(mean_num)])
disp(['Mean dist. =  ' num2str(mean_d)])

%%
figure
for i = 1:numel(dist)
    plot(dist(i).num,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
title('$Num.(d<5)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/frame$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$N_{tracers}$','interpreter','latex',FontWeight='bold',FontSize=18)
plot(xlim,[mean_num mean_num],'k--',LineWidth=2)
text(0,mean_num*4,'$\overline{N_{tracer}} = 1.4$','interpreter','latex',FontWeight='bold',FontSize=18)

grid on
axis padded

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*pNexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.');
set(col,'position',[0.68 0.88 0.2 0.02])

% V
figure
tf=tiledlayout(3,1,'TileSpacing','compact');
nexttile
for i = 1:numel(dist)
    plot(dist(i).Vrelx,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
ylabel('$V^{rel}_x(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*pNexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.');
set(col,'position',[0.68 0.88 0.2 0.02])

nexttile
for i = 1:numel(dist)
    plot(dist(i).Vrely,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
ylabel('$V^{rel}_y(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on

nexttile
for i = 1:numel(dist)
    plot(dist(i).Vrelz,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
ylabel('$V^{rel}_z(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/frame$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on

% A
figure
tf=tiledlayout(3,1,'TileSpacing','compact');
nexttile
for i = 1:numel(dist)
    plot(dist(i).Vrelx,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
ylabel('$A^{rel}_x(mm/s^2)$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*pNexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.');
set(col,'position',[0.68 0.88 0.2 0.02])

nexttile
for i = 1:numel(dist)
    plot(dist(i).Arely,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
ylabel('$A^{rel}_y(mm/s^2)$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on

nexttile
for i = 1:numel(dist)
    plot(dist(i).Arelz,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
ylabel('$A^{rel}_z(mm/s^2)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/frame$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on