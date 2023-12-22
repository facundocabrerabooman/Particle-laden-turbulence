function plotLagr_pdfs_lavision(pdfV,pdfA,mycolormap,ifsave)

color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
lw= 1;

figure;

semilogy(pdfV(1).xpdf,pdfV(1).pdf,'d-',MarkerSize=5,Color=color3(1,:),LineWidth=lw);hold on;
semilogy(pdfV(2).xpdf,pdfV(2).pdf,'d-',MarkerSize=5,Color=color3(2,:),LineWidth=lw);
semilogy(pdfV(3).xpdf,pdfV(3).pdf,'d-',MarkerSize=5,Color=color3(3,:),LineWidth=lw);

xpdfn.V1 = linspace(-1.2,1.0,1024)*1e3;
xpdfn.V2 = linspace(-1.0,1.0,1024)*1e3;
xpdfn.V3 = linspace(-1.1,0.8,1024)*1e3;
semilogy(xpdfn.V1,normpdf(xpdfn.V1,pdfV(1).mean,pdfV(1).std),'--',Color=color3(1,:),LineWidth=lw);
semilogy(xpdfn.V2,normpdf(xpdfn.V2,pdfV(2).mean,pdfV(2).std),'--',Color=color3(2,:),LineWidth=lw);
semilogy(xpdfn.V3,normpdf(xpdfn.V3,pdfV(3).mean,pdfV(3).std),'--',Color=color3(3,:),LineWidth=lw);

% lavision output
% m/s to mm/s
load('LavisionOutput\tracers.mat')
[pdfVL.x,xpdfVL.x] = hist(d(:,5),256); 
[pdfVL.y,xpdfVL.y] = hist(d(:,6),256); 
[pdfVL.z,xpdfVL.z] = hist(d(:,7),256); 
semilogy(xpdfVL.x.*1e3,pdfVL.x/sum(pdfVL.x),'-',Color=color3(1,:),LineWidth=lw); hold on
semilogy(xpdfVL.y.*1e3,pdfVL.y/sum(pdfVL.y),'-',Color=color3(2,:),LineWidth=lw); 
semilogy(xpdfVL.z.*1e3,pdfVL.z/sum(pdfVL.z),'-',Color=color3(3,:),LineWidth=lw); 

set(gca,FontSize=15)
legend('$V_x$','$V_y$','$V_z$','$Gfitx$','$Gfity$','$Gfitz$','$V_xL$','$V_yL$','$V_zL$','interpreter','latex',Location='best',FontSize=12);
title('$PDFn$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(V)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$V(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

if ifsave==1
    savefig('./Lagrangian Figures/PDFnsV')
    saveas(gcf,'./Lagrangian Figures/PDFnsV','png')
    saveas(gcf,'./Lagrangian Figures/PDFnsV','pdf')
end

%%
figure
semilogy(pdfA(1).xpdf,pdfA(1).pdf,'^-',MarkerSize=5,Color=color3(1,:),LineWidth=lw);hold on
semilogy(pdfA(2).xpdf,pdfA(2).pdf,'^-',MarkerSize=5,Color=color3(2,:),LineWidth=lw);
semilogy(pdfA(3).xpdf,pdfA(3).pdf,'^-',MarkerSize=5,Color=color3(3,:),LineWidth=lw);

xpdfn.A1 = linspace(-0.8e5,0.8e5,1024);
xpdfn.A2 = linspace(-0.8e5,0.8e5,1024);
xpdfn.A3 = linspace(-0.8e5,0.8e5,1024);
semilogy(xpdfn.A1,normpdf(xpdfn.A1,pdfA(1).mean,pdfA(1).std),'--',MarkerSize=5,Color=color3(1,:),LineWidth=lw);
semilogy(xpdfn.A2,normpdf(xpdfn.A2,pdfA(2).mean,pdfA(2).std),'--',MarkerSize=5,Color=color3(2,:),LineWidth=lw);
semilogy(xpdfn.A3,normpdf(xpdfn.A3,pdfA(3).mean,pdfA(3).std),'--',MarkerSize=5,Color=color3(3,:),LineWidth=lw);


% lavision output
% m/s to mm/s
[pdfAL.x,xpdfAL.x] = hist(d(:,8),256); 
[pdfAL.y,xpdfAL.y] = hist(d(:,9),256); 
[pdfAL.z,xpdfAL.z] = hist(d(:,10),256); 
semilogy(xpdfAL.x*1e3,pdfAL.x/sum(pdfAL.x),'-',Color=color3(1,:),LineWidth=lw); hold on
semilogy(xpdfAL.y*1e3,pdfAL.y/sum(pdfAL.y),'-',Color=color3(2,:),LineWidth=lw); 
semilogy(xpdfAL.z*1e3,pdfAL.z/sum(pdfAL.z),'-',Color=color3(3,:),LineWidth=lw); 

set(gca,FontSize=15)
legend('$A_x$','$A_y$','$A_z$','$Gfitx$','$Gfity$','$Gfitz$','$A_xL$','$A_yL$','$A_zL$','interpreter','latex',Location='best',FontSize=12);
title('$PDFn$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(A)n$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$A(mm/s^2)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded
% xlim([-5 5])

if ifsave==1
    savefig('./Lagrangian Figures/PDFnsA')
    saveas(gcf,'./Lagrangian Figures/PDFnsA','png')
    saveas(gcf,'./Lagrangian Figures/PDFnsA','pdf')
end