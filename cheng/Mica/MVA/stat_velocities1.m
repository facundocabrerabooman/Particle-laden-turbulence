function [stat]=stat_velocities(vel,n)

%c=340;
%fact=c/2/nu0/sin(theta/2)*32768;
%fact=1;
[B,A]=butter(4,10*2/32768);
data=[];
acc=[];
velf=[];
nbulle=[];
% ll=[];
% ll=sum(vel.length(vel.good))
% data=zeros(1,ll);
% jj=1;
for j=1:length(vel.good)
    
    %if((vel.status(j)==1)&(max(abs(vel.data(vel.good(j)).seg).^2)>1e-6))
        %data(jj:jj+vel.length(vel.good(j))-1)=vel.data(vel.good(j)).freq*fact;
        data=[data vel.data(vel.good(j)).freq];
        %data=[data vel.data(vel.good(j)).velf*fact];
        nbulle=[nbulle ones(1,length(vel.data(vel.good(j)).freq))*j];
        if isfield(vel.data,'acc');
            acc=[acc vel.data(vel.good(j)).acc];
        end
        if isfield(vel.data,'velf');
            velf=[velf vel.data(vel.good(j)).velf];
        end
    %end
    %jj=jj+vel.length(vel.good(j));
end

%figure;plot(data);
stat.N=numel(data);
stat.Nbulles=numel(vel.good);
stat.mean=mean(data);
stat.std=std(data);
stat.flatness=kurtosis(data);
stat.skewness=skewness(data);
[stat.hist stat.x]=hist(data,64);
[stat.xpdf,stat.pdf]=mkpdf3(data,64,-5,5);
if isfield(vel.data,'acc');
    stat.acc.N=numel(acc);
    stat.acc.mean=mean(acc);
    stat.acc.std=std(acc);
    
    std(acc-stat.acc.mean)/stat.acc.std
    
    stat.acc.flatness=kurtosis(acc);
    stat.acc.skewness=skewness(acc);
    [x,pdf]=mkpdf3(acc,64,-10,10);
    stat.acc.xpdf=x;
    stat.acc.pdf=pdf;
end
if isfield(vel.data,'velf');
    stat.velf.N=numel(velf);
    stat.velf.mean=mean(velf);
    stat.velf.std=std(velf);
    
    std(velf-stat.velf.mean)/stat.velf.std
    
    stat.velf.flatness=kurtosis(velf);
    stat.velf.skewness=skewness(velf);
    [x,pdf]=mkpdf3(velf,64,-10,10);
    stat.velf.xpdf=x;
    stat.velf.pdf=pdf;
end
