function sig=sig_sup(vel_in,varargin)

% sigs=sig_sup(vel_in)
%
% calcul sigma_sup(tau), la variance de vitesse des segments plus longs que
% tau. le signal de sortie a N+1 elements (l'element 1 correspondant a
% tau=0 et l'element N+1 a tau=N)
%

% load(vel_in)
% eval(['vel=' vel_in '.data(' vel_in '.good);'])
% eval(['clear ' vel_in]);

vel=vel_in;
clear vel_in;

if nargin==1
    L=arrayfun(@(x)(numel(x.velf)),vel.data(vel.good));
    Nmin=min(L);
    Nmax=max(L);
    mean_tot=sum(arrayfun(@(x)(mean(x.velf)),vel.data(vel.good)).*L)/sum(L);
    sig_tot=sum(arrayfun(@(x)(mean((x.velf-mean_tot).^2)),vel.data(vel.good)).*L)/sum(L);
    
    sig(1:Nmin+1)=sig_tot;
    for k=Nmin+1:Nmax
        I=find(L>=k);
        sig(k+1)=sum(arrayfun(@(x)(mean((x.velf-mean_tot).^2)),vel.data(vel.good(I))).*L(I))/sum(L(I));
    end
else
    vfield=varargin{1};
    L=arrayfun(@(x)(numel(x.(vfield))),vel);
    Nmin=min(L);
    Nmax=max(L);
    mean_tot=0.*sum(arrayfun(@(x)(mean(x.(vfield))),vel).*L)/sum(L);
    sig_tot=sum(arrayfun(@(x)(mean((x.(vfield)-mean_tot).^2)),vel).*L)/sum(L);
    
    sig(1:Nmin+1)=sig_tot;
    for k=Nmin+1:Nmax
        I=find(L>=k);
        sig(k+1)=sum(arrayfun(@(x)(mean((x.(vfield)-mean_tot).^2)),vel(I)).*L(I))/sum(L(I));
    end
    sig(end)=[];
end
