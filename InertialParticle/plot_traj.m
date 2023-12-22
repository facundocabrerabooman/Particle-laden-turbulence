function plot_traj(part,titlestring,axisrange,ifsave,fout)

%%
fields = fieldnames(part);
for i = 1:length(fields)
    part2.(fields{i}) = vertcat(part.(fields{i}));
end

if ~exist(fout)
    mkdir(fout)
end

scale = 0.008;
%%
% axisrange = [-5 25 -18 10 -40 40];

Nframemax = 1e6;
nexp= fix(part2.T/Nframemax)+1;
uni_exp = unique(nexp);
Nexp = numel(uni_exp);

mycolormap = mycolor('#B10000','#FFFFFF','#0000B2');
mycolormap2 = mycolor('#063970','#eeeee4','#e28743');
mycolorind = 1:floor(size(mycolormap,1)/Nexp):size(mycolormap,1);
%% all trajectories 3D plot
figure;
scatter3(part2.Xf,part2.Zf,part2.Yf,abs(part2.Ay)*scale,part2.Vy,'filled',LineWidth=2)
axis equal;box;grid on;axis padded
set(gca,FontSize=15)
legend('$size \propto Acce. $','interpreter','latex',FontWeight='bold',FontSize=15)
title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
zlabel('$y(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$z(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$x(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
legend(Location='northeast')
colormap(mycolormap);col =colorbar;
ylabel(col,'$Vy(m/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
axis(axisrange)

if ifsave == 1
    savefig([fout '\traj3D'])
    saveas(gcf,[fout,'\traj3D'],'png')
    saveas(gcf,[fout,'\traj3D'],'pdf')
end
%% all trajectories 2D plot - xz projection 
figure;
scatter(part2.Xf,part2.Zf,abs(part2.Ay)*scale,part2.Vy,'filled',LineWidth=2);hold on
scatter(part2.Xf(1),part2.Zf(1),abs(part2.Ay(1))*scale*0.2,'go',MarkerFaceColor='g',LineWidth=2)
axis equal;box;grid;axis padded
set(gca,FontSize=15)
legend('$size \propto Acce. $','interpreter','latex',FontWeight='bold',FontSize=15)
title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$z(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$x(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap);col =colorbar;
ylabel(col,'$Vy(m/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
axis(axisrange(1:4))

if ifsave ==1
    savefig([fout '\traj2Dxz'])
    saveas(gcf,[fout,'\traj2Dxz'],'png')
    saveas(gcf,[fout,'\traj2Dxz'],'pdf')
end
%% xz projection
figure;
for i = 1:Nexp
    plot(part2.Xf(nexp==i),part2.Zf(nexp==i),Color=mycolormap2(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end
plot(part2.Xf(1),part2.Zf(1),'go',MarkerFaceColor='g',LineWidth=2);hold on
axis equal;box on;grid on;axis padded
set(gca,FontSize=15)
title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$z(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$x(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',ceil([0.01 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.')
set(col,'position',[0.65 0.2 0.2 0.02])
axis(axisrange(1:4))

if ifsave == 1
    savefig([fout '\traj2Dxz2'])
    saveas(gcf,[fout,'\traj2Dxz2'],'png')
    saveas(gcf,[fout,'\traj2Dxz2'],'pdf')
end