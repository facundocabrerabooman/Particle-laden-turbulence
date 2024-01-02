%extraction de frequence par methode MVA avec filtre de kalman
%
% parametres
%
% NbMic taille de la fenetre temporelle
% K nombre de fenetre successives utilises pour constuire la matrice de correlation
% NbS nombre de sources sonores recherchees
% Teta vecteur colonne des valeurs initiales
% MaxIter nb max d'iterations effedctues pour converger
% epsilon valeur du test de convergence
% Reps est la matrice de covariance du bruit de modelisation
%
% les frequences sont dans la matrice freq
% iHes est une estimation de la qualite de l'estimation
% Puis est la matrice des amplitudes des sources

%MVAinit=0.9;

NbMic=7;
K=13;
NbS=1;
Teta=[2*pi*MVAinit];  %valeurs d'initialisation
MaxIter=20;  %nb max d'iterations
epsilon=1e-8; %test de convergence
Reps=1e-5*eye(NbS);

[l,ll]=size(sig);
if ll>l
    sig=sig.';
end;

tot=length(sig)-K-NbMic-1;
freq=zeros(NbS,tot);
iHes=zeros(NbS,tot);
syg=zeros(1,tot);
Rx=zeros(NbMic,NbMic);
vp=zeros(NbMic,tot);
Puis=zeros(NbS,tot);




%**********************************************************************
%                PREMIERE ESTIMATION
%**********************************************************************

ii=1;

% Acquisition d'une nouvelle matrice Rx
Rx=matcor(sig(ii:ii+K+NbMic-2),NbMic);

% Decomposition en elts propres, calcul de PiB et Ry.   

[U,VP,V]=svd(Rx);
VP=diag(VP);
vp(:,ii)=VP;

PiB=zeros(NbMic,NbMic);
syg2=0;
for k=(NbS+1):NbMic;
    PiB=PiB+U(:,k)*U(:,k)';
    syg2=syg2+VP(k);
end;
syg2=abs(syg2/(NbMic-NbS));		% valeur estimee du bruit
Ry=zeros(NbMic,NbMic);
for k=1:NbS;
    Ry=Ry+(VP(k)-syg2)*U(:,k)*U(:,k)';
end; 

Nb_Iter=0;
Veps=10;

while ((Nb_Iter<MaxIter)&((Veps'*Veps)>epsilon^2));     % condition de CV
    
    Nb_Iter=Nb_Iter+1;																						
    Tetac=Teta;                              
    
    Sc=zeros(NbMic,NbS);
    for k=1:NbS;
        Sc(:,k)=exp(1i*Tetac(k)*(0:(NbMic-1))')/sqrt(NbMic);
    end;                                     
    DSc=diag(1i*(-(NbMic-1)/2:1:(NbMic-1)/2))*Sc;  
    iSSt=inv(Sc'*Sc);
    PiBc=eye(NbMic)-Sc*iSSt*Sc';
    P=iSSt*Sc'*Ry*Sc*iSSt;
    A=DSc'*PiBc*PiB;
    Grad=(2*K/syg2)*real(diag(A*Sc*P));
    Hessian=(2*K/syg2)*real((A*PiBc*DSc).*conj(P));
    iHess=inv(Hessian);
    Veps=iHess*Grad;                        
    Teta=Tetac-Veps';
end;

syg(ii)=sqrt(syg2);
iHes(:,ii)=diag(iHess);
freq(:,ii)=Teta.'/2/pi;
Puis(1:NbS,ii)=diag(P);


%*****************************************************************
%                   ESTIMATIONS SUIVANTES
%*****************************************************************

Hess0=Hessian;
Gamma=iHess;

progress=0;
for ii=2:tot
    %if ii-fix(ii/100)*100==0 disp((100*ii/tot));end;
    if 100*ii/tot-progress > 10 
        disp((100*ii/tot));
        progress=100*ii/tot;
    end;
    
    % Acquisition d'une nouvelle matrice Rx
    Rx=matcor(sig(ii:ii+K+NbMic-2),NbMic);
    
    % Decomposition en elts propres, calcul de PiB et Ry.   
    [U,VP,V]=svd(Rx);
    VP=diag(VP);
    vp(:,ii)=VP;
    
    PiB=zeros(NbMic,NbMic);
    syg2=0;
    for k=(NbS+1):NbMic;
        PiB=PiB+U(:,k)*U(:,k)';
        syg2=syg2+VP(k);
    end;
    syg2=abs(syg2/(NbMic-NbS));		% valeur estimee du bruit
    Ry=zeros(NbMic,NbMic);
    for k=1:NbS;
        Ry=Ry+(VP(k)-syg2)*U(:,k)*U(:,k)';
    end; 
    %Ry=Rx; 
    
 
    Tetac=Teta;                              
    
    Sc=zeros(NbMic,NbS);
    for k=1:NbS
        Sc(:,k)=exp(1i*Tetac(k)*(0:(NbMic-1))');%/sqrt(NbMic);
    end;                                     
    DSc=diag(1i*(-(NbMic-1)/2:1:(NbMic-1)/2))*Sc;  
    iSSt=inv(Sc'*Sc);
    PiBc=eye(NbMic)-Sc*iSSt*Sc';
    P=iSSt*Sc'*Ry*Sc*iSSt;
    A=DSc'*PiBc*PiB;
    Grad=(2*K/syg2)*real(diag(A*Sc*P));
    Hessian=(2*K/syg2)*real((A*PiBc*DSc).*conj(P));
    iHess=inv(Hessian);
    %Veps=iHess*Grad;
    iHes(:,ii)=diag(iHess);                    % mise � jour du hessien
    Gammac=Gamma+Reps;
    Hess0=inv(Gammac)+Hessian;
    iHess=inv(Hess0);
    Gamma=iHess;
    
    % nouvelle estimation de teta
    Teta=Tetac-(iHess*Grad)';
    
    syg(ii)=sqrt(syg2);
    freq(:,ii)=Teta.'/2/pi;

    Puis(1:NbS,ii)=diag(P);
end;

freq=[freq;iHes;real(Puis)];


