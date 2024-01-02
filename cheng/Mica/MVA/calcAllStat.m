function stat=calcAllstat(directory)
dd=dir([directory '/*G1']);
dd(1)='';
dd(1)='';
N=numel(dd);
Npdf=129;
for k=1:numel(dd)
	if dd(k).isdir==1
		cd(dd(k).name);
		param=load('param');
		dv=dir('stat*');
		S=load(dv.name);
		F=fieldnames(S);
		
		dv=dir('velm*');
		vel=load(dv.name);
		Fvel=fieldnames(vel);
		vel=accelerations(vel.(char(Fvel)),30);
		S.(char(F))=stat_velocities5(vel,Npdf);
		%S.(char(F))=stat_velocities5(vel.(char(Fvel)),Npdf);
		
		eval([char(F) '=S.(char(F));']);
		eval(['save ' char(F) ' ' char(F) ';']);
		
		%stat(k).param=param;
		%stat(k).dirname=dd(k).name;
		%stat(k).stat=S.(char(F));
		
		cd ..;
	end
end

stat=mergeStat(directory);