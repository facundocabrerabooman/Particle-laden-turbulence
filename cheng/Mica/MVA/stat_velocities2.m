function [stat]=stat_velocities2(vel,n,varargin)

% stat=stat_velocities2(vel,n)
%
% calcul des stattistiques de la serie de mesures de vitesse Lagrangienne
% de la structure vel.
% - calcul des pdfs normalisees et non normalisees de vitesse, vitesse filtree 
%   et acceleration
% - estimation a partir des pdf des moment d'ordre 1,2,3 et 4

if nargin >2
    Xv=varargin{1};
    Xvf=varargin{2};
    Xa=varargin{3};
    if numel(Xv)~=n
        disp(sprintf('verifier la taille de Xv'));
        return
    end
    if numel(Xv)~=n
        disp(sprintf('verifier la taille de Xvf'));
        return
    end
    if numel(Xv)~=n
        disp(sprintf('verifier la taille de Xa'));
        return
    end
else
    wv=15;
    wvf=15;
    wa=15;

    data=[];
    acc=[];
    velf=[];
    nbulle=[];
    % ll=[];
    % ll=sum(vel.length(vel.good))
    % data=zeros(1,ll);
    % jj=1;

    for j=1:min(100,numel(vel.good))
        data=[data vel.data(vel.good(j)).velf];
        if isfield(vel.data,'acc');
            acc=[acc vel.data(vel.good(j)).acc];
        end
        if isfield(vel.data,'velf');
            velf=[velf vel.data(vel.good(j)).velf];
        end
    end
    stat.mean=mean(data);
    stat.std=std(data);
    %stat.flatness=kurtosis(data);
    %stat.skewness=skewness(data);
    %[stat.hist stat.x]=hist(data,64);
    %[stat.xpdf,stat.pdf]=mkpdf3(data,64,-5,5);
    if isfield(vel.data,'acc');

        stat.acc.mean=mean(acc);
        stat.acc.std=std(acc);

        %   std(acc-stat.acc.mean)/stat.acc.std

        %   stat.acc.flatness=kurtosis(acc);
        %   stat.acc.skewness=skewness(acc);
        %   [x,pdf]=mkpdf3(acc,64,-10,10);
        %   stat.acc.xpdf=x;
        %   stat.acc.pdf=pdf;
    end
    if isfield(vel.data,'velf');

        stat.velf.mean=mean(velf);
        stat.velf.std=std(velf);

        %  std(velf-stat.velf.mean)/stat.velf.std

        %  stat.velf.flatness=kurtosis(velf);
        %  stat.velf.skewness=skewness(velf);
        %  [x,pdf]=mkpdf3(velf,64,-10,10);
        %  stat.velf.xpdf=x;
        %  stat.velf.pdf=pdf;
    end

    Xv=linspace(stat.mean-wv*stat.std,stat.mean+wv*stat.std,n);
    Xvf=linspace(stat.velf.mean-wvf*stat.velf.std,stat.velf.mean+wvf*stat.velf.std,n);
    Xa=linspace(stat.acc.mean-wa*stat.acc.std,stat.acc.mean+wa*stat.acc.std,n);
end

Nv=zeros(1,n);
Nvf=zeros(1,n);
Na=zeros(1,n);
Lv=0;
Lvf=0;
La=0;

for j=1:length(vel.good)
        %if((vel.status(j)==1)&(max(abs(vel.data(vel.good(j)).seg).^2)>1e-6))
        %data(jj:jj+vel.length(vel.good(j))-1)=vel.data(vel.good(j)).freq*fact;
        %data=[data vel.data(vel.good(j)).freq*fact];
        %data=[data vel.data(vel.good(j)).velf*fact];
        %nbulle=[nbulle ones(1,length(vel.data(vel.good(j)).freq))*j];
        
        Nv=Nv+hist(vel.data(vel.good(j)).freq,Xv);
        Lv=Lv+numel(vel.data(vel.good(j)).freq);
        
        if isfield(vel.data,'acc');
            Na=Na+hist(vel.data(vel.good(j)).acc,Xa);
            La=La+numel(vel.data(vel.good(j)).acc);
            %acc=[acc vel.data(vel.good(j)).acc*fact];
        end
        if isfield(vel.data,'velf');
            Nvf=Nvf+hist(vel.data(vel.good(j)).velf,Xvf);
            Lvf=Lvf+numel(vel.data(vel.good(j)).velf);
            %velf=[velf vel.data(vel.good(j)).velf*fact];
        end
    %end
    %jj=jj+vel.length(vel.good(j));
end
%stat.std

stat.N=Lv;
stat.velf.N=Lvf;
stat.acc.N=La;
stat.Nbulles=numel(vel.good);

stat.pdf=Nv/integ(Nv,Xv);
stat.xpdf=Xv;
stat.mean=integ(stat.xpdf.*stat.pdf,stat.xpdf);
stat.std=sqrt(integ((stat.xpdf-stat.mean).^2.*stat.pdf,stat.xpdf));
stat.skewness=integ((stat.xpdf-stat.mean).^3.*stat.pdf,stat.xpdf)/stat.std^3;
stat.flatness=integ((stat.xpdf-stat.mean).^4.*stat.pdf,stat.xpdf)/stat.std^4;

if isfield(vel.data,'velf');
    stat.velf.pdf=Nvf/integ(Nvf,Xvf);
    stat.velf.xpdf=Xvf;
    stat.velf.mean=integ(stat.velf.xpdf.*stat.velf.pdf,stat.velf.xpdf);
    stat.velf.std=sqrt(integ((stat.velf.xpdf-stat.velf.mean).^2.*stat.velf.pdf,stat.velf.xpdf));
    stat.velf.skewness=integ((stat.velf.xpdf-stat.velf.mean).^3.*stat.velf.pdf,stat.velf.xpdf)/stat.velf.std^3;
    stat.velf.flatness=integ((stat.velf.xpdf-stat.velf.mean).^4.*stat.velf.pdf,stat.velf.xpdf)/stat.velf.std^4;
end
if isfield(vel.data,'acc');
    stat.acc.pdf=Na/integ(Na,Xa);
    stat.acc.xpdf=Xa;
    stat.acc.mean=integ(stat.acc.xpdf.*stat.acc.pdf,stat.acc.xpdf);
    stat.acc.std=sqrt(integ((stat.acc.xpdf-stat.acc.mean).^2.*stat.acc.pdf,stat.acc.xpdf));
    stat.acc.skewness=integ((stat.acc.xpdf-stat.acc.mean).^3.*stat.acc.pdf,stat.acc.xpdf)/stat.acc.std^3;
    stat.acc.flatness=integ((stat.acc.xpdf-stat.acc.mean).^4.*stat.acc.pdf,stat.acc.xpdf)/stat.acc.std^4;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDF  Normalisees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nv=zeros(1,n);
Nvf=zeros(1,n);
Na=zeros(1,n);
if nargin>2
    Xv=(Xv-stat.mean)/stat.std;
    Xvf=(Xvf-stat.velf.mean)/stat.velf.std;
    Xa=(Xa-stat.acc.mean)/stat.acc.std;
else
    Xv=linspace(-wv,wv,n);
    Xvf=linspace(-wvf,wvf,n);
    Xa=linspace(-wa,wa,n);
end

for j=1:length(vel.good)
        Nv=Nv+hist((vel.data(vel.good(j)).freq-stat.mean)/stat.std,Xv);
        if isfield(vel.data,'acc');
            Na=Na+hist((vel.data(vel.good(j)).acc-stat.acc.mean)/stat.acc.std,Xa);
        end
        if isfield(vel.data,'velf');
            Nvf=Nvf+hist((vel.data(vel.good(j)).velf-stat.velf.mean)/stat.velf.std,Xvf);
        end
end
stat.pdfn=Nv/integ(Nv,Xv);
stat.xpdfn=Xv;

if isfield(vel.data,'velf');
    stat.velf.pdfn=Nvf/integ(Nvf,Xvf);
    stat.velf.xpdfn=Xvf;
end
if isfield(vel.data,'acc');
    stat.acc.pdfn=Na/integ(Na,Xa);
    stat.acc.xpdfn=Xa;
end