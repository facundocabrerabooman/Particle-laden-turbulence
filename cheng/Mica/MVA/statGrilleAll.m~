function sGrille=statGrilleAll(rep,prefixe)


c=340;
nu0=80e3;
theta=165*pi/180;

xm=2.5;
xM=3.25;

X=linspace(xm,xM,32);
x0=(xm+xM)/2;

fact=c/2/nu0/sin(theta/2)*32768;

fstat=find_params([],rep,[prefixe '*.mat']);

N=0;
for k=1:numel(fstat)
	Nd=find(fstat{k}==filesep);
	D=fstat{k}(1:(Nd(numel(Nd))-1));
	Ns=strfind(fstat{k},[filesep prefixe]);
	suffixe=fstat{k}(Ns+numel(prefixe)+1:numel(fstat{k}));
	
	var=who('-FILE',fstat{k},'stat*');
	load(fstat{k},'stat*');

	stat_tmp=eval(char(var));
	clear(char(var),'var');
	
	if isfield(stat_tmp,'vel')
		U=stat_tmp.vel.mean*fact;
	elseif isfield(stat_tmp,'velf');
		U=stat_tmp.velf.mean*fact;
	else
		U=inf;
	end
	
	if U~=inf
		N=N+1;
		grille=paramGrille(X,U);
		ff=fieldnames(grille);
		mgrille=structfun(@(x)(mean(x)),grille);
		dgrille=structfun(@(x)(std(x)),grille);
		mgrille=cell2struct(num2cell(mgrille),ff);
		dgrille=cell2struct(num2cell(dgrille),ff);

		statGrille.statm_file=fstat{k};
		statGrille.mgrille=mgrille;
		statGrille.dgrille=dgrille;

		sGrille(N)=statGrille;
		
		fgrille=[D filesep prefixe '_Grille' suffixe];
		save(fgrille,'statGrille');

	end
end
	


% --------------------------------------------------------------------    
function params=find_params(params,rep,ext)
param_file=dir(fullfile(rep,ext));
if ~isempty(param_file)
	
    for kk=1:numel(param_file)
        params_tmp=[rep filesep param_file(kk).name];
        params=[params ; cellstr(params_tmp)];
    end
end

dd=dir(rep);
reps=[];

for kk=1:numel(dd)
    if(dd(kk).isdir==1)&((strcmp(dd(kk).name,'.')==0)&(strcmp(dd(kk).name,'..')==0))
        cd(dd(kk).name);     
        params=find_params(params,pwd,ext);
        cd(rep);
    end
end