function stat=mergeStat(directory)
dd=dir(directory);
dd(1)='';
dd(1)='';
N=numel(dd);

for k=1:numel(dd)
	if dd(k).isdir==1
		cd(dd(k).name);
		param=load('param');
		dv=dir('stat*');
		S=load(dv.name);
		F=fieldnames(S);
		stat(k).param=param;
		stat(k).dirname=dd(k).name;
		stat(k).stat=S.(char(F));
		
		cd ..;
	end
end