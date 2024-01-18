function Sn=structFunc_struct(S,field,n)

disp('Structure function for each trajectory...');
[sn,weights]=arrayfun(@(x)(structFunc(x.(field),n)),S,'UniformOutput',false);

disp('Global structure function...');
sn=cell2struct(sn,'Sn');
weights=cell2struct(weights,'w');
Sn=meanstruct(sn,'Sn',weights);
Sn.tau=0:(numel(Sn.mean)-1);

end