function [x,pdf]=mkpdf3(s,taille,a,b,varargin)

% [x,pdf]=mkpdf3(s,taille,a,b)

if nargin>4
    norm_flag=varargin{1};
else
    norm_flag=0;
end


if norm_flag==1
    sn=(s-mean(s))/std(s);
elseif norm_flag==2
    sn=s/std(s);
elseif norm_flag==3
    sn=s-mean(s);
else
    sn=s;
end

[n x]=hist(sn,linspace(a,b,taille));
pdf=n/(sum(n)*(x(2)-x(1)));
i1=findi(x,a);
i2=findi(x,b);
x=x(i1:i2);
pdf=pdf(i1:i2);

% pdf=medfilt1(pdf,8);
% figure(1);
% subplot(1,2,1)
% plot(x,pdf);
% hold on
% %plot(x0,pdf0,'r');
% grid
% hold off
% subplot(1,2,2)
% plot(x,log(pdf))
% hold on
% %plot(x0,log(pdf0),'r');
% grid
% hold off

