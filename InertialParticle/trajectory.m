clear all;clc
close all

fin = 'F:\MaltlabOutput\0811_1g_noTurb_4kfps\particle';
fout = [fin '\Figures'];

%%
load([fin '\STB.mat'])

%%
titlestring = '$1g-without Turb.$';
axisrange = [-5 25 -18 10 -40 40];

fps = 4230;
Nframemax = 1e6;
NTrackmax = 1e9;
nexp= fix(part.T/Nframemax)+1;
t = mod(part.T,Nframemax)/fps;
uni_exp = unique(nexp);
Nexp = numel(uni_exp);

mycolormap1 = mycolor('#B10000','#FFFFFF','#0000B2');
mycolormap2 = mycolor('#063970','#eeeee4','#e28743');
mycolorind = 1:floor(size(mycolormap1,1)/Nexp):size(mycolormap1,1);
%% all trajectories 3D plot
figure;
scatter3(part.X,part.Z,part.Y,abs(part.Ay)*0.5,part.Vy,'filled',LineWidth=2)
axis equal;box;grid on;axis padded
set(gca,FontSize=15)
legend('$size \propto Acce. $','interpreter','latex',FontWeight='bold',FontSize=15)
title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
zlabel('$y(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$z(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$x(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
legend(Location='northeast')
colormap(mycolormap1);col =colorbar;
ylabel(col,'$Vy(m/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
axis(axisrange)

savefig([fout '\traj3D'])
saveas(gcf,[fout,'\traj3D'],'png')
saveas(gcf,[fout,'\traj3D'],'pdf')

%% all trajectories 2D plot - xz projection 
figure;
scatter(part.X,part.Z,abs(part.Ay)*0.5,part.Vy,'filled',LineWidth=2);hold on
scatter(part.X(1),part.Z(1),abs(part.Ay(1))*0.5,'g*',LineWidth=2)
axis equal;box;grid;axis padded
set(gca,FontSize=15)
legend('$size \propto Acce. $','interpreter','latex',FontWeight='bold',FontSize=15)
title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$z(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$x(mm)$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap1);col =colorbar;
ylabel(col,'$Vy(m/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
axis(axisrange(1:4))

savefig([fout '\traj2Dxz'])
saveas(gcf,[fout,'\traj2Dxz'],'png')
saveas(gcf,[fout,'\traj2Dxz'],'pdf')

%% xz projection
figure;
for i = 1:Nexp
    plot(part.X(nexp==i),part.Z(nexp==i),Color=mycolormap2(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end
plot(part.X(1),part.Z(1),'g*',LineWidth=2);hold on
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

savefig([fout '\traj2Dxz2'])
saveas(gcf,[fout,'\traj2Dxz2'],'png')
saveas(gcf,[fout,'\traj2Dxz2'],'pdf')
