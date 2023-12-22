function [stat]=stat_velocities3(vel,n,varargin)

% stat=stat_velocities3(vel,n,varargin)
%
% varargin{1}= 1 --> moyenne des pdfs ponderee par nbre de points par serie
%            = 0 --> meme poids pour chaque serie independamment du nbre de points  
%
% calcul des statistiques de la serie de mesures de vitesse Lagrangienne
% de la structure vel.
% - calcul des pdfs normalisees et non normalisees de vitesse, vitesse filtree 
%   et acceleration
% - estimation a partir des pdf des moment d'ordre 1,2,3 et 4
%
% la difference avec stat_velocities2 reside dans le fait que la version 3
% calcule les pdfs pour chaque serie de mesure independamment puis estime
% la pdf globale comme la moyenne des pdfs, alors que la version 2 calcule
% une pdf a partir de toutes les realisations de toutes les series.

if nargin>2
    flag_pond=varargin{1};
else
    flag_pond=0;
end

vel_s=split_vel(vel);

stat.xpdft=zeros(1,n);
stat.pdft=zeros(1,n);
stat.pdfnt=zeros(1,n);;
stat.xpond=zeros(numel(vel_s),n);

stat.velf.xpdft=zeros(1,n);
stat.velf.pdft=zeros(1,n);
stat.velf.pdfnt=zeros(1,n);
stat.velf.xpond=zeros(numel(vel_s),n);

stat.acc.xpdft=zeros(1,n);
stat.acc.pdft=zeros(1,n);
stat.acc.pdfnt=zeros(1,n);
stat.acc.xpond=zeros(numel(vel_s),n);

for jj=1:numel(vel_s);
    if nargin>3
        stat_tmp=stat_velocities2(vel_s(jj),n,varargin{2:4});
    else
        if jj==1    
            stat_tmp=stat_velocities2(vel_s(jj),n);
        else
            stat_tmp=stat_velocities2(vel_s(jj),n,stat.pdf(1).xpdf,stat.velf.pdf(1).xpdf,stat.acc.pdf(1).xpdf);
        end
    end
    
    stat.N(jj)=stat_tmp.N;
    stat.Nbulles(jj)=stat_tmp.Nbulles;
    stat.mean(jj)=stat_tmp.mean;
    stat.std(jj)=stat_tmp.std;
    stat.skewness(jj)=stat_tmp.skewness;
    stat.flatness(jj)=stat_tmp.flatness;
    stat.pdf(jj).xpdf=stat_tmp.xpdf;
    stat.pdf(jj).pdf=stat_tmp.pdf;
    stat.pdf(jj).xpdfn=stat_tmp.xpdfn;
    stat.pdf(jj).pdfn=stat_tmp.pdfn;
    if flag_pond>0
        stat.pond(jj)=stat.N(jj);
    else
        stat.pond(jj)=1;
    end
    if flag_pond==2
        stat.xpond(jj,:)=sign(stat.pdf(jj).pdf)*stat.pond(jj);
    else
        stat.xpond(jj,:)=ones(1,n)*stat.pond(jj);
    end
    
    stat.pdft=stat.pdft+stat.pdf(jj).pdf*stat.pond(jj);
    stat.pdfnt=stat.pdft+stat.pdf(jj).pdfn*stat.pond(jj);
    
    if isfield(vel.data,'velf');
        stat.velf.N(jj)=stat_tmp.velf.N;
        stat.velf.mean(jj)=stat_tmp.velf.mean;
        stat.velf.std(jj)=stat_tmp.velf.std;
        stat.velf.skewness(jj)=stat_tmp.velf.skewness;
        stat.velf.flatness(jj)=stat_tmp.velf.flatness;
        stat.velf.pdf(jj).xpdf=stat_tmp.velf.xpdf;
        stat.velf.pdf(jj).pdf=stat_tmp.velf.pdf;
        stat.velf.pdf(jj).xpdfn=stat_tmp.velf.xpdfn;
        stat.velf.pdf(jj).pdfn=stat_tmp.velf.pdfn;
        if flag_pond>0
            stat.velf.pond(jj)=stat.velf.N(jj);
        else
            stat.velf.pond(jj)=1;
        end
        if flag_pond==2
            stat.velf.xpond(jj,:)=sign(stat.velf.pdf(jj).pdf)*stat.velf.pond(jj);
        else
            stat.velf.xpond(jj,:)=ones(1,n)*stat.velf.pond(jj);
        end
        stat.velf.pdft=stat.velf.pdft+stat.velf.pdf(jj).pdf*stat.velf.pond(jj);
        stat.velf.pdfnt=stat.velf.pdfnt+stat.velf.pdf(jj).pdfn*stat.velf.pond(jj);
    end
    if isfield(vel.data,'acc');
        stat.acc.N(jj)=stat_tmp.N;
        stat.acc.mean(jj)=stat_tmp.acc.mean;
        stat.acc.std(jj)=stat_tmp.acc.std;
        stat.acc.skewness(jj)=stat_tmp.acc.skewness;
        stat.acc.flatness(jj)=stat_tmp.acc.flatness;
        stat.acc.pdf(jj).xpdf=stat_tmp.acc.xpdf;
        stat.acc.pdf(jj).pdf=stat_tmp.acc.pdf;
        stat.acc.pdf(jj).xpdfn=stat_tmp.acc.xpdfn;
        stat.acc.pdf(jj).pdfn=stat_tmp.acc.pdfn;
        if flag_pond>0
            stat.acc.pond(jj)=stat.acc.N(jj);
        else
            stat.acc.pond(jj)=1;
        end
        
        if flag_pond==2
            stat.acc.xpond(jj,:)=sign(stat.acc.pdf(jj).pdf)*stat.acc.pond(jj);
        else
            stat.acc.xpond(jj,:)=ones(1,n)*stat.acc.pond(jj);
        end
        stat.acc.xpond(jj,:)=sign(stat.acc.pdf(jj).pdf)*stat.acc.pond(jj);
        stat.acc.pdft=stat.acc.pdft+stat.acc.pdf(jj).pdf*stat.acc.pond(jj);
        stat.acc.pdfnt=stat.acc.pdfnt+stat.acc.pdf(jj).pdfn*stat.acc.pond(jj);
    end
    
end

stat.xpdft=stat.pdf(1).xpdf;
stat.xpdfnt=stat.pdf(1).xpdfn;
N=sum(stat.xpond,1);
[y,i]=find(N~=0);
stat.pdft(i)=stat.pdft(i)./N(i);
stat.pdfnt(i)=stat.pdfnt(i)./N(i);

if isfield(vel.data,'velf');
    stat.velf.xpdft=stat.velf.pdf(1).xpdf;
    stat.velf.xpdfnt=stat.velf.pdf(1).xpdfn;
    N=sum(stat.velf.xpond,1);
    [y,i]=find(N~=0);
    stat.velf.pdft(i)=stat.velf.pdft(i)./N(i);
    stat.velf.pdfnt(i)=stat.velf.pdfnt(i)./N(i);
end

if isfield(vel.data,'acc');
    stat.acc.xpdft=stat.acc.pdf(1).xpdf;
    stat.acc.xpdfnt=stat.acc.pdf(1).xpdfn;
    N=sum(stat.acc.xpond,1);
    [y,i]=find(N~=0);
    stat.acc.pdft(i)=stat.acc.pdft(i)./N(i);
    stat.acc.pdfnt(i)=stat.acc.pdfnt(i)./N(i);
end