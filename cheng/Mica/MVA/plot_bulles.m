function [bulles]=plot_bulles(vel,i1,N)

hold on
data=[];
[B,A]=butter(4,.005);

xt=vel.length(i1)/2;
for j=i1:i1+N
        data=[data filtfilt(B,A,abs(vel.data(j).seg)'.^2)];
        %data=[data abs(vel.data(j).seg)'.^2];
        %data=[data medfilt1(abs(vel.data(j).seg)'.^2,400)];
        text(xt,max(abs(vel.data(j).seg)'.^2)*1.2,num2str(j));
        xt=xt+vel.length(j)/2+vel.length(j+1)/2;

end
plot(data,'r');
hold off;