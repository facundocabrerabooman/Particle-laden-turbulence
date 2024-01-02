function [Rub,Sub]=autocorr_struct_unbiased(S1,field1)

%
% [Rub,Sub]=autocorr_struct_unbiased(S1,field1)
%
% Calcule la correlation Rub et la fonction de structure Sub debiasees de
% telle sorte que la relation Rub=1-Sub/2sigma soit verifiée.
%


% First, estimate the standard correlation and structure functions...
ll=arrayfun(@(x)(numel(x.(field1))),S1);
moy=sum(arrayfun(@(x)(mean(x.(field1))),S1).*ll)/sum(ll);
std2=sum(arrayfun(@(x)(mean((x.(field1)-moy).^2)),S1).*ll)/sum(ll)


disp(sprintf('Moyennes bulle a bulle ...\n'));
[r, weights_ub]=arrayfun(@(x)(autocorr(x.(field1),moy)),S1,'UniformOutput',false);
[sn, weights_ub]=arrayfun(@(x)(structFunc(x.(field1),2)),S1,'UniformOutput',false);

disp(sprintf('Moyenne d''ensemble ...\n'));
r=cell2struct(r,'R');
sn=cell2struct(sn,'Sn');
weights_ub=cell2struct(weights_ub,'w');

R=meanstruct(r,'R',weights_ub);
Sn=meanstruct(sn,'Sn',weights_ub);

% Then, unbiase ...
fc=(2*R.mean+Sn.mean)/2/R.mean(1);
Rub=1-(1-R.mean/R.mean(1))./fc;
Sub=2*R.mean(1)*(1-fc)./fc+Sn.mean./fc;
