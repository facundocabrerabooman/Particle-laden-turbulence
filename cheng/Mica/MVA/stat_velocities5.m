function [stat]=stat_velocities5(velm,n,varargin)

% stat=stat_velocities5(vel,n)
%
% calcul des statistiques de la serie de mesures de vitesse Lagrangienne
% de la structure vel.
% - calcul des pdfs normalisees et non normalisees de vitesse, vitesse filtree 
%   et acceleration
% - estimation a partir des pdf des moment d'ordre 1,2,3 et 4
% - estimation directe des moments d'ordre 1,2,3 et 4 (par arrayfun)


L=arrayfun(@(x)(numel(x.velf)),velm.data(velm.good));


vel.mean=sum(arrayfun(@(x)(mean(x.velf)),velm.data(velm.good)).*L)/sum(L);
vel.std=sqrt(sum(arrayfun(@(x)(mean((x.velf-vel.mean).^2)),velm.data(velm.good)).*L)/sum(L));
vel.meanstd=sum(arrayfun(@(x)(std(x.velf)),velm.data(velm.good))).*L)/sum(L);
vel.tauIndiv=sum(arrayfun(@(x)(std(x.velf)./mean(x.velf)),velm.data(velm.good)).*L)/sum(L);
vel.skewness=sum(arrayfun(@(x)(mean((x.velf-vel.mean).^3)),velm.data(velm.good)).*L)/sum(L)/vel.std^3;
vel.flatness=sum(arrayfun(@(x)(mean((x.velf-vel.mean).^4)),velm.data(velm.good)).*L)/sum(L)/vel.std^4;
vel.pdf=mkpdf5(velm.data(velm.good),'velf',n);

acc.mean=sum(arrayfun(@(x)(mean(x.acc)),velm.data(velm.good)).*L)/sum(L);
acc.std=sqrt(sum(arrayfun(@(x)(mean((x.acc-acc.mean).^2)),velm.data(velm.good)).*L)/sum(L));
acc.meanstd=sum(arrayfun(@(x)(std(x.acc)),velm.data(velm.good)).*L)/sum(L);
acc.tauIndiv=sum(arrayfun(@(x)(std(x.acc)./mean(x.acc)),velm.data(velm.good)).*L)/sum(L);
acc.skewness=sum(arrayfun(@(x)(mean((x.acc-acc.mean).^3)),velm.data(velm.good)).*L)/sum(L)/acc.std^3;
acc.flatness=sum(arrayfun(@(x)(mean((x.acc-acc.mean).^4)),velm.data(velm.good)).*L)/sum(L)/acc.std^4;
acc.pdf=mkpdf5(velm.data(velm.good),'acc',n);

%Nrj=cell2struct(arrayfun(@(x)((x.velf-vel.mean).*(x.acc-acc.mean)),velm.data(velm.good)(velm.good),'UniformOutput',false),'nrj')';
Nrj=cell2struct(arrayfun(@(x)((x.velf).^2),velm.data(velm.good)(velm.good),'UniformOutput',false),'nrj')';

nrj.mean=sum(arrayfun(@(x)(mean(x.nrj)),Nrj).*L)/sum(L);
nrj.std=sqrt(sum(arrayfun(@(x)(mean((x.nrj-nrj.mean).^2)),Nrj).*L)/sum(L));
nrj.meanstd=sum(arrayfun(@(x)(std(x.nrj)),Nrj).*L)/sum(L);
nrj.tauIndiv=sum(arrayfun(@(x)(std(x.nrj)./mean(x.nrj)),Nrj).*L)/sum(L);
nrj.skewness=sum(arrayfun(@(x)(mean((x.nrj-nrj.mean).^3)),Nrj).*L)/sum(L)/nrj.std^3;
nrj.flatness=sum(arrayfun(@(x)(mean((x.nrj-nrj.mean).^4)),Nrj).*L)/sum(L)/nrj.std^4;
nrj.pdf=mkpdf5(Nrj,'nrj',n);

Nrj2=cell2struct(arrayfun(@(x)((x.velf-vel.mean).*(x.acc-acc.mean)),velm.data(velm.good)(velm.good),'UniformOutput',false),'nrj')';

nrj2.mean=sum(arrayfun(@(x)(mean(x.nrj)),Nrj2).*L)/sum(L);
nrj2.std=sqrt(sum(arrayfun(@(x)(mean((x.nrj-nrj2.mean).^2)),Nrj2).*L)/sum(L));
nrj2.meanstd=sum(arrayfun(@(x)(std(x.nrj)),Nrj2).*L)/sum(L);
nrj2.tauIndiv=sum(arrayfun(@(x)(std(x.nrj)./mean(x.nrj)),Nrj2).*L)/sum(L);
nrj2.skewness=sum(arrayfun(@(x)(mean((x.nrj-nrj2.mean).^3)),Nrj2).*L)/sum(L)/nrj2.std^3;
nrj2.flatness=sum(arrayfun(@(x)(mean((x.nrj-nrj2.mean).^4)),Nrj2).*L)/sum(L)/nrj2.std^4;
nrj2.pdf=mkpdf5(Nrj,'nrj',n);


stat.vel=vel;
stat.acc=acc;
stat.nrj=nrj;
stat.nrj2=nrj2;


