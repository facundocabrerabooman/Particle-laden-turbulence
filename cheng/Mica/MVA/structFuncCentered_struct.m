function Sn=structFuncCentered_struct(S1,field1,n)
%
% Sn=SructFunc(S1,field1,n)
% 
% calculates Sn centered with forward and backwart time lag



disp(sprintf('Moyennes bulle a bulle ...\n'));
[snbwd, snfwd, weights_ub]=arrayfun(@(x)(structFuncCentered(x.(field1),n)),S1,'UniformOutput',false);

disp(sprintf('Moyenne d''ensemble ...\n'));
sn=cell2struct(snbwd,'Sn');
weights_ub=cell2struct(weights_ub,'w');
SnBwd=meanstruct(sn,'Sn',weights_ub);
SnBwd.tau=-[0:(numel(SnBwd.mean)-1)];


sn=cell2struct(snfwd,'Sn');
SnFwd=meanstruct(sn,'Sn',weights_ub);
SnFwd.tau=0:(numel(SnBwd.mean)-1);

Sn.Bwd = SnBwd;
Sn.Fwd = SnFwd;