function plotLagr_pdfs(pdfV,pdfA,mycolormap,ifsave,fout)

color1 = '#a155b9';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw1= 2;
lw2 = 1;

figure;
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=lw1);hold on;
semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=lw1);
semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=lw1);

semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'o-',MarkerSize=5,Color=color3(1,:),LineWidth=lw2);
semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'o-',MarkerSize=5,Color=color3(2,:),LineWidth=lw2);
semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'o-',MarkerSize=5,Color=color3(3,:),LineWidth=lw2);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=lw2);

set(gca,FontSize=15)
legend('$V_x$','$V_y$','$V_z$','$A_x$','$A_y$','$A_z$','interpreter','latex',Location='best',FontSize=12);
title('$PDF$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(V,A)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$V, A$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% add subfigure
% axes('Position',[0.22 0.62 0.22 0.22]);
% semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',MarkerSize=2,Color=color3(1,:),LineWidth=lw);hold on;
% semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',MarkerSize=2,Color=color3(2,:),LineWidth=lw);
% semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',MarkerSize=2,Color=color3(3,:),LineWidth=lw);
% 
% semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'^-',MarkerSize=2,Color=color3(1,:),LineWidth=lw);
% semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'^-',MarkerSize=2,Color=color3(2,:),LineWidth=lw);
% semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'^-',MarkerSize=2,Color=color3(3,:),LineWidth=lw);
% 
% xpdf=linspace(-5,5,1024);
% plot(xpdf,normpdf(xpdf,0,1),'k',LineWidth=lw);
% grid on
% set(gca,FontSize=12)
% xlim([-5 5])

if ifsave ==1
    figname = ['./' fout '/PDFs'];
    savefig(figname)
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'pdf')
%     save('LagrangianStat.mat','pdfV','pdfA','-append')
end
