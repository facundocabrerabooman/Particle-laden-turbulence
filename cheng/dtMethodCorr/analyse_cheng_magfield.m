load tracks.mat
% keep only tracks longer than 10 frames
var=tracklong;
fps=3000;% sampling freq
clear track Ilong Ine tracklong
% var = 
%   1×45797 struct array with fields:
%     T
%     Ntrack
%     X
%     Y
%     Z
%     Xf
%     Yf
%     Zf
%     Vx
%     Vy
%     Vz
%     Ax
%     Ay
%     Az
%     Tf
%     Ntrackf

%%
ltrack=zeros(1,length(var));
for k=1:length(var)
    ltrack(k)=length(var(k).T);
end

%%
h=figure;
axesh = axes('Parent',h);
plot(ltrack,'bo');grid
set(axesh,'FontSize',16);
clear axesh h
%%

nmax=10;
clear ndt md* meand* count
ndt=1:nmax;
mdx=zeros(length(var),nmax);
mdx2=zeros(length(var),nmax);
md2x=zeros(length(var),nmax);
md2x2=zeros(length(var),nmax);
count=zeros(1,nmax);


for k=1:length(var)
    x=var(k).X;
        for kn=1:nmax;
        n=kn;% number of frames used in computation of dx and d2x
            if length(x)>2*n && sum(isnan(x))==0
            dx=dx_over_n_points(x,n);
            d2x=d2x_over_n_points(x,n);

            count(kn)=count(kn)+1;
            mdx(k,kn)=mean(dx);
            mdx2(k,kn)=mean(dx.^2);
            md2x(k,kn)=mean(d2x);
            md2x2(k,kn)=mean(d2x.^2);
            end
            clear dx d2x;
        end
    clear x
    if round(k/5000)*5000==k;
    disp(k);
    end
end

meandx=sum(mdx,1)./count;
meandx2=sum(mdx2,1)./count;
meand2x=sum(md2x,1)./count;
meand2x2=sum(md2x2,1)./count;

clear md*
%%
h=figure;
axesh = axes('Parent',h);
plot(ndt,meandx,'bo');grid
set(axesh,'FontSize',16);
clear axesh h
xlabel('$ndt$','interpreter','latex','VerticalAlignment','top');
ylabel('$\langle dx \rangle$','interpreter','latex','VerticalAlignment','bottom');
title('mean displacement')

h=figure;
axesh = axes('Parent',h);
plot(ndt.^2,meandx2,'bo');grid
set(axesh,'FontSize',16);
clear axesh h
xlabel('$ndt^2$','interpreter','latex','VerticalAlignment','top');
ylabel('$\langle dx^2 \rangle$','interpreter','latex','VerticalAlignment','bottom');
title('mean squared displacement')

h=figure;
axesh = axes('Parent',h);
plot(ndt.^4,meand2x2,'bo');grid
set(axesh,'FontSize',16);
clear axesh h
xlabel('$ndt^4$','interpreter','latex','VerticalAlignment','top');
ylabel('$\langle d^2x^2 \rangle$','interpreter','latex','VerticalAlignment','bottom');
title('mean squared second order increment')
% compute increments by taking average 
