function find_optimal_ndt(track,nmaxdt,nmaxtau,fieldin,ifsave,vargin)
%
% find_optimal_ndt
% show plots of <dx(t)dx(t+tau)> and <d2x(t)d2x(t+tau)>
% so that to select the optimal starting points and the length of the ndt
% used for calculating the correlations
% 
% inputs: 
% tracks: contains the fields 'X','Y, 'Z' (raw positions)
% nmaxdt: max number of ndt
% nmaxtau: max number of ndtau
% fieldin: 'X','Y, 'Z';
% ifsavefig: save flag
% 
% Cheng Wang 21/09/23



% clear;clc
% load tracks.mat
% % keep only tracks longer than 10 frames
% var=tracklong;
% fps=3000;% sampling freq
% clear long_thres wopt lopt tracklong
% addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));
%%
% nmaxdt=10;
% nmaxtau = 10;
ndt=1:nmaxdt;
ntau = 1:nmaxtau;

mdxdtau  = zeros(nmaxtau,nmaxdt,length(track));
md2xdtau = zeros(nmaxtau,nmaxdt,length(track));
counttau = zeros(nmaxtau,nmaxdt);

for k=1:length(track)
    x=track(k).(fieldin);
        for kdt=1:nmaxdt
            n=kdt;% number of frames used in computation of dx and d2x
            if length(x)>2*n
                dx=dx_over_n_points(x,n);
                d2x=d2x_over_n_points(x,n);

                for ktau = 1:nmaxtau
                    nn = ktau;
                    if nn<(length(dx)-1) && nn<(length(d2x)-1)
                        counttau(ktau,kdt)=counttau(ktau,kdt)+1;
                        mdxdtau(ktau,kdt,k) = mean(dx(1:end-nn).*dx(1+nn:end));
                        md2xdtau(ktau,kdt,k) =mean(d2x(1:end-nn).*d2x(1+nn:end));
                    end
                end
            end
            clear dx d2x;
        end
    clear x
end
%%
meandxdtau2 = sum(mdxdtau,3)./counttau;
meand2xdtau2 = sum(md2xdtau,3)./counttau;

%% find optimal ndt for fixed ntau
mycolormap = mycolor('#063970','#eeeee4','#e28743');
mycolorind1 = round(linspace(1,size(mycolormap,1),nmaxtau));
mycolorind2 = round(linspace(1,size(mycolormap,1),nmaxdt));

f1 = figure;
tiledlayout(2,1)
nexttile
for i = 1:nmaxtau
    plot(ndt.^2,meandxdtau2(i,:),'-o',Color=mycolormap(mycolorind1(i),:),LineWidth=2);hold on
end
grid
xlabel('$(n \cdot dt)^2$','interpreter','latex','VerticalAlignment','top');
% xticks(ndt.^2)
% xticklabels(string(ndt))
ylabel('$\langle dx(t) dx(t+\tau) \rangle$','interpreter','latex','VerticalAlignment','bottom');
title('Choose your Optimal ndt')
set(gca,'FontSize',16);

nexttile
for i = 1:nmaxtau
    plot(ndt.^4,meand2xdtau2(i,:),'-o',Color=mycolormap(mycolorind1(i),:),LineWidth=2);hold on
end
grid
xlabel('$(n \cdot dt)^4$','interpreter','latex','VerticalAlignment','top');
ylabel('$\langle d^2x(t) d^2x(t+\tau) \rangle$','interpreter','latex','VerticalAlignment','bottom');
set(gca,'FontSize',16);
colormap(mycolormap);
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',ceil([0.01 0.25 0.5 0.75 1]*nmaxtau),FontSize=12,Location='eastoutside');
xlabel(col,'$n \cdot \tau$','interpreter','latex',FontWeight='bold',FontSize=18)

%% check dependence on tau for different ndt
f2 = figure;
tiledlayout(2,1)
nexttile
for i = 1:nmaxdt
    plot(ntau,meandxdtau2(:,i)/meandxdtau2(1,i),'-o',Color=mycolormap(mycolorind2(i),:),LineWidth=2);hold on
%     plot(ntau,meandxdtau2(:,i),'-o',Color=mycolormap(mycolorind2(i),:),LineWidth=2);hold on
end
grid
xlabel('$n \cdot \tau$','interpreter','latex','VerticalAlignment','top');
ylabel('$\frac{\langle dx(t) dx(t+\tau)}{(dx(t))^2 \rangle}$','interpreter','latex','VerticalAlignment','bottom');
title('Check dependence on tau for different ndt')
set(gca,'FontSize',16);

nexttile
for i = 1:nmaxdt
    plot(ntau,meand2xdtau2(:,i)/meand2xdtau2(1,i),'-+',Color=mycolormap(mycolorind2(i),:),LineWidth=2);hold on
%     plot(ntau,meand2xdtau2(:,i),'-+',Color=mycolormap(mycolorind2(i),:),LineWidth=2);hold on
end
grid
xlabel('$n \cdot \tau$','interpreter','latex','VerticalAlignment','top');
ylabel('$\frac{\langle d^2x(t) d^2x(t+\tau)}{(d^2x(t))^2 \rangle}$','interpreter','latex','VerticalAlignment','bottom');
set(gca,'FontSize',16);
colormap(mycolormap);
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',ceil([0.01 0.25 0.5 0.75 1]*nmaxdt),FontSize=12,Location='eastoutside');
xlabel(col,'$n \cdot dt$','interpreter','latex',FontWeight='bold',FontSize=18)

%%
if strcmp(ifsave,'save') == 1
    fout = [vargin 'Figures_dtCorr'];
    if ~exist(fout) == 1
        mkdir(fout)
        disp(['mkdir ... ' fout])
    end
    figname1 = [fout '/dtCorr_dt_' fieldin];
    savefig(f1, figname1)
    saveas(f1,figname1,'png')
    saveas(f1,figname1,'pdf')
    figname2 = [fout '/dtCorr_tau_' fieldin];
    savefig(f2, figname2)
    saveas(f2,figname2,'png')
    saveas(f2,figname2,'pdf')
end
%%
% meandxdtau = sum(mdxdtau,[2,3])./sum(counttau,2);
% meand2xdtau = sum(md2xdtau,[2,3])./sum(counttau,2);
% 
% figure;plot(ntau,meandxdtau,'bo');grid
% xlabel('$n\tau$','interpreter','latex','VerticalAlignment','top');
% ylabel('$\langle dx(t) dx(t+\tau) \rangle$','interpreter','latex','VerticalAlignment','bottom');
% title('mean displacement - tau')
% set(gca,'FontSize',16);
% 
% figure;plot(ntau,meand2xdtau,'bo');grid
% xlabel('$n\tau$','interpreter','latex','VerticalAlignment','top');
% ylabel('$\langle d^2x(t) d^2x(t+\tau) \rangle$','interpreter','latex','VerticalAlignment','bottom');
% title('mean squared displacement - tau')
% set(gca,'FontSize',16);

