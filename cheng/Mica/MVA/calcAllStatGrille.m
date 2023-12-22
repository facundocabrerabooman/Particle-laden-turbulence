function stat=calcAllStatGrille(directory)

fact=70.232851783117198;
dd=dir(directory);
dd(1)='';
dd(1)='';
N=numel(dd);
Npdf=129;
h=waitbar(0,'Calcul des parametres de la grille ...');

for k=1:numel(dd)
	waitbar(k/numel(dd),h);
	if dd(k).isdir==1
		cd(dd(k).name);
		param=load('param');
		dv=dir('stat*');
		S=load(dv.name);
		F=fieldnames(S);
		
		%dv=dir('velm*');
		%vel=load(dv.name);
		%Fvel=fieldnames(vel);
		
		%S.(char(F))=stat_velocities5(vel.(char(Fvel)),Npdf);
		
		paramGrille=calcParamGrille(linspace(2.5,3.25,64),S.(char(F)).vel.mean*fact);
		

		paramGrille.depsilon=std(paramGrille.epsilon);
		paramGrille.epsilon=mean(paramGrille.epsilon);

		paramGrille.deta=std(paramGrille.eta);
		paramGrille.eta=mean(paramGrille.eta);
		
		paramGrille.dtau_eta=std(paramGrille.tau_eta);
		paramGrille.tau_eta=mean(paramGrille.tau_eta);
		
		paramGrille.dtau=std(paramGrille.tau);
		paramGrille.tau=mean(paramGrille.tau);
		
		paramGrille.dL=std(paramGrille.L);		
		paramGrille.L=mean(paramGrille.L);
		
		paramGrille.dRlambda=std(paramGrille.Rlambda);		
		paramGrille.Rlambda=mean(paramGrille.Rlambda);
		
		paramGrille.dlambda=std(paramGrille.lambda);	
		paramGrille.lambda=mean(paramGrille.lambda);
		
		paramGrille.dTe=std(paramGrille.Te);
		paramGrille.Te=mean(paramGrille.Te);
		
		S.(char(F)).paramGrille=paramGrille;
		
		eval([char(F) '=S.(char(F));']);
		eval(['save ' char(F) ' ' char(F) ';']);
		
		cd ..;
	end
end
close(h);

stat=mergeStat(directory);