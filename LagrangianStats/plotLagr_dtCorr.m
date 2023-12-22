function plotLagr_dtCorr(tau,corrv,corra,ifsave,fout)

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];

fields = fieldnames(tau);
figure
subplot(2,1,1)
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau.(f),corrv.(f)/corrv.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end
legend('$X$','$Y$','$Z$','interpreter','latex',Location='best',FontSize=12);
xlabel('$\tau(s)$','interpreter','latex',FontSize=18);
ylabel('$\langle u(t)u(t+\tau) \rangle$','interpreter','latex',FontSize=18);
grid;set(gca,'FontSize',15);

subplot(2,1,2)
for kfield = 1:numel(fields)
    f = fields{kfield};
    plot(tau.(f),corra.(f)/corra.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
end
% plot(xlim,[0 0])
grid;set(gca,'FontSize',15);
legend('$X$','$Y$','$Z$','interpreter','latex',Location='best',FontSize=12);
xlabel('$\tau(s)$','interpreter','latex',FontSize=18);
ylabel('$\langle a(t)a(t+\tau) \rangle$','interpreter','latex',FontSize=18);


if ifsave ==1
    figname = ['./' fout '/dtCorr'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
end