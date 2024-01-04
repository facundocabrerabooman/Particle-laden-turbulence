function [tau,corrv,corra] = dtCorr(track,varargin)

% dtCorr(track) or dtCorr(track,ndts,ndtl,fps)
% returns the denoised correlation functions of
% derivates (velocity, acceleration) using dt method
% 
% inputs: 
% tracks: contains the fields 'X','Y, 'Z' (raw positions)
% ndts: start point for ndt, ndts = [5 5 5];
% ndtl: length for ndt, ndtl = [10 10 10];
% fps: frame rate, fps = 2996;
% 
% outputs:
% tau
% corrv: correalation functions of velocity
% corra: correalation functions of acceleartions
%
% Cheng Wang 21/09/23

if isempty(varargin)
    ndts = [6 6 6];
    ndtl = [10 10 10];
    fps = 2996;
else
    ndts = varargin{1};
    ndtl = varargin{2};
    fps = varargin{3};
end

ndt.X = ndts(1):1:(ndts(1)+ndtl(1)-1);
ndt.Y = ndts(2):1:(ndts(2)+ndtl(2)-1);
ndt.Z = ndts(2):1:(ndts(3)+ndtl(3)-1);

Lmax = 2000;
Ndtmax = max(ndtl);
Nmaxtrack = numel(track);

fields = {'X','Y','Z'};
%% mean
for kfield = 1:numel(fields)
    f = fields{kfield};
    for ktrack = 1:Nmaxtrack
        t = track(ktrack).(f);
        for kdt = 1:numel(ndt.(f))
            n = ndt.(f)(kdt);
            if length(t)>2*n
                dx = dx_over_n_points(t,n);
                d2x = d2x_over_n_points(t,n);
                
                L.(f)(ktrack,kdt) = length(dx);
                % mean of each track
                mdx.(f)(ktrack,kdt) = mean(dx);
                md2x.(f)(ktrack,kdt) = mean(d2x);
                clear dx d2x
            end
        end
    end
    % global mean
    mean_mdx.(f) = sum(mdx.(f).*L.(f),1)./sum(L.(f),1);
    mean_md2x.(f) = sum(md2x.(f).*L.(f),1)./sum(L.(f),1);
end


%% substact mean
for kfield = 1:numel(fields)
    f = fields{kfield}
    mdx0 = mean_mdx.(f);
    md2x0 = mean_md2x.(f);
    cx=zeros(2*Lmax-1,Ndtmax);
    c2x=zeros(2*Lmax-1,Ndtmax);
    nc=zeros(2*Lmax-1,Ndtmax);
    for ktrack = 1:Nmaxtrack
        t = track(ktrack).(f);
        for kdt = 1:numel(ndt.(f))
            n = ndt.(f)(kdt);
            if length(t)>2*n
                dx = dx_over_n_points(t,n);
                d2x = d2x_over_n_points(t,n);
                
                % substact mean value
                dx = dx-mdx0(kdt);
                d2x = d2x-md2x0(kdt);
                ll = length(dx);
                
                cx(Lmax-ll+1:Lmax+ll-1,kdt) = cx(Lmax-ll+1:Lmax+ll-1,kdt) + xcorr(dx,'none');
                c2x(Lmax-ll+1:Lmax+ll-1,kdt) = c2x(Lmax-ll+1:Lmax+ll-1,kdt) + xcorr(d2x,'none');
                nc(Lmax-ll+1:Lmax+ll-1,kdt) = nc(Lmax-ll+1:Lmax+ll-1,kdt) + [1:ll,ll-1:-1:1]';
    
                clear dx d2x
            end
        end
    end
    mcx.(f) = cx./nc;
    mc2x.(f) = c2x./nc;

%     figure
%     for i =1:Ndtmax
%         idx = ~isnan(mcx.(f)(:,i));
%         plot(mcx.(f)(idx,i));hold on
%     end
%     figure
%     for i =1:Ndtmax
%         idx = ~isnan(mc2x.(f)(:,i));
%         plot(mc2x.(f)(idx,i));hold on
%     end

    tau0.(f)=(-Lmax+1:Lmax-1)'/fps;
    
    kmax = size(mcx.(f),1);
    corrv0.(f) = zeros(kmax,2);
    corra0.(f) = zeros(kmax,2);
    for i = 1:kmax
        corrv0.(f)(i,:)= polyfit(ndt.(f).^2,mcx.(f)(i,:),1);
        corra0.(f)(i,:)= polyfit(ndt.(f).^4,mc2x.(f)(i,:),1);
    end
    corrv0.(f)(:,2) = [];
    corra0.(f)(:,2) = [];
    
    idx = find(~isnan(corrv0.(f)));
    tau1.(f) = tau0.(f)(idx);
    corrv1.(f) = corrv0.(f)(idx);
    corra1.(f) = corra0.(f)(idx);
    
    tau.(f)=tau1.(f)((length(tau1.(f))-1)/2+1:end);
    corrv.(f)=corrv1.(f)((length(corrv1.(f))-1)/2+1:end);
    corra.(f)=corra1.(f)((length(corra1.(f))-1)/2+1:end);
end
%%
% mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
% color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
% 
% fields = fieldnames(tau);
% figure
% subplot(2,1,1)
% for kfield = 1:numel(fields)
%     f = fields{kfield};
%     plot(tau.(f),corrv.(f)/corrv.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
% end
% legend('$X$','$Y$','$Z$','interpreter','latex',Location='best',FontSize=12);
% xlabel('$\tau(s)$','interpreter','latex',FontSize=18);
% ylabel('$\langle u(t)u(t+\tau) \rangle$','interpreter','latex',FontSize=18);
% grid;set(gca,'FontSize',15);
% 
% subplot(2,1,2)
% for kfield = 1:numel(fields)
%     f = fields{kfield};
%     plot(tau.(f),corra.(f)/corra.(f)(1),'-',Color=color3(kfield,:),LineWidth = 2);hold on
% end
% plot(xlim,[0 0])
% grid;set(gca,'FontSize',15);
% legend('$X$','$Y$','$Z$','interpreter','latex',Location='best',FontSize=12);
% xlabel('$\tau(s)$','interpreter','latex',FontSize=18);
% ylabel('$\langle a(t)a(t+\tau) \rangle$','interpreter','latex',FontSize=18);
% 
% figname = ['./' fout '/dtCorr'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')