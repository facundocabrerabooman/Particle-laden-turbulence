clear
clc
close all

%%
fps = 4230;

ppart = load('.\particle\STB.mat');
tpart = load('.\tracers\STB.mat');

%%
Nframemax = 1e6;
NTrackmax = 1e9;

emax = 150;
dist(1:emax) = struct('numIn5',[],'dIn5',[],'numIn10',[],'dIn10',[]);
%%
pnexp = fix(ppart.part.T/Nframemax)+1;
pt = mod(ppart.part.T,Nframemax)/fps;
puni_exp = unique(pnexp);
pNexp = numel(puni_exp);

%%
tnexp = fix(tpart.part.T/Nframemax)+1;
tt = mod(tpart.part.T,Nframemax)/fps;
tuni_exp = unique(tnexp);
tNexp = numel(tuni_exp);

%%
mycolormap = mycolor('#063970','#eeeee4','#e28743');%('#063970','#eeeee4','#e28743')
mycolorind = 1:floor(size(mycolormap,1)/pNexp):size(mycolormap,1);
patchcolor1 = mycolormap(size(mycolormap,1)/2,:);
patchcolor2 = mycolor('#a155b9');
patchalpha1 = 0.4;
patchalpha2 = 0.3;

%%
for i = 1:pNexp
    pind = find(pnexp == puni_exp(i));
    ptime = pt(pind);

    tind = find(tnexp == tuni_exp(i));
    ttime = tt(tind);

    for j = 1:numel(ptime)
        k=find(ttime==ptime(j));
        d{j,:} = sqrt((ppart.part.X(pind(j)) - tpart.part.X(tind(k))).^2 + ...
                    (ppart.part.X(pind(j)) - tpart.part.X(tind(k))).^2 + ...
                    (ppart.part.X(pind(j)) - tpart.part.X(tind(k))).^2);
    
        ind1 = find(d{j,:}<5);
        dist(i).numIn5(j) = size(ind1,1);
        dist(i).dIn5{j} = d{j,:}(ind1);

        ind2 = find(d{j,:}<10);
        dist(i).numIn10(j) = size(ind2,1);
        dist(i).dIn10{j}= d{j,:}(ind2);
    end
    
end
dist = dist(arrayfun(@(X)(~isempty(X(1).numIn5~=0)),dist));


%%
mean_numIn5  = mean([dist.numIn5]);
mean_dIn5  = mean(arrayfun(@(X)(mean(vertcat(X.dIn5{1,:}))),dist));
disp(['Mean num in 5mm neighbor  =  ' num2str(mean_numIn5)])
disp(['Mean dist. in 5mm neighbor  =  ' num2str(mean_dIn5)])

mean_numIn10 = mean([dist.numIn10]);
mean_dIn10 = mean(arrayfun(@(X)(mean(vertcat(X.dIn10{1,:}))),dist));
disp(['Mean num in 10mm neighbor = ' num2str(mean_numIn10)])
disp(['Mean dist. in 10mm neighbor = ' num2str(mean_dIn10)])

%%
figure
tf=tiledlayout(2,1,'TileSpacing','compact');
nexttile
for i = 1:numel(dist)
    plot(dist(i).numIn5,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
title('$Num.(d<5)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/frame$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$N_{tracers}$','interpreter','latex',FontWeight='bold',FontSize=18)
plot(xlim,[mean_numIn5 mean_numIn5],'k--',LineWidth=2)
text(0,mean_numIn5*0.3,'$\overline{N_{tracer}} = 38.2$','interpreter','latex',FontWeight='bold',FontSize=18)

grid on
axis padded

nexttile
for i = 1:numel(dist)
    for j = 1:1
% for i = 1:1
%     for j = 1:numel(dist(i).dIn5)
        plot(dist(i).dIn5{1,j},'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
    end
end

set(gca,FontSize=15)
title('$In\ the \ first \ frame:\ |X_{t}-X{p}|(d<5) $','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$N_{tracers}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$d/mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
plot(xlim,[mean_dIn5 mean_dIn5],'k--',LineWidth=2)
text(0,mean_dIn5*0.3,'$\overline{d} = 2.44$','interpreter','latex',FontWeight='bold',FontSize=18)

grid on
axis padded

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*pNexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.');
set(col,'position',[0.68 0.88 0.2 0.02])

savefig('./Num5')
saveas(gcf,'./Num5','png')
saveas(gcf,'./Num5','pdf')

%%
figure
tf=tiledlayout(2,1,'TileSpacing','compact');
nexttile
for i = 1:numel(dist)
    plot(dist(i).numIn10,'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
end
set(gca,FontSize=15)
title('$Num.(d<10)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/frame$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$N_{tracers}$','interpreter','latex',FontWeight='bold',FontSize=18)
plot(xlim,[mean_numIn10 mean_numIn10],'k--',LineWidth=2)
text(0,mean_numIn10*0.4,'$\overline{N_{tracer}} = 70.2$','interpreter','latex',FontWeight='bold',FontSize=18)

grid on
axis padded

nexttile
for i = 1:numel(dist)
    for j = 1:1
% for i = 1:1
%     for j = 1:numel(dist(i).dIn5)
        plot(dist(i).dIn10{1,j},'-',Color=mycolormap(mycolorind(i),:),LineWidth=2);hold on
    end
end

set(gca,FontSize=15)
title('$In\ the \ first \ frame:\ |X_{t}-X{p}|(d<10) $','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$N_{tracers}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$d/mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
plot(xlim,[mean_dIn10 mean_dIn10],'k--',LineWidth=2)
text(0,mean_dIn10*0.3,'$\overline{d} = 4.71$','interpreter','latex',FontWeight='bold',FontSize=18)

grid on
axis padded

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*pNexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.');
set(col,'position',[0.68 0.88 0.2 0.02])

savefig('./Num10')
saveas(gcf,'./Num10','png')
saveas(gcf,'./Num10','pdf')

