function [Sn Sig Sigdt]=xcorr_struct(S1,field1,n)

% Rxx=xcorr_struct(S1,field1,n)
%



disp(sprintf('Moyennes bulle a bulle ...\n'));
[sn, weights_ub, sigma, sigmadt]=arrayfun(@(x)(xcorrFunc(x.(field1),n)),S1,'UniformOutput',false);

disp(sprintf('Moyenne d''ensemble ...\n'));
sn=cell2struct(sn,'Sn');
weights_ub=cell2struct(weights_ub,'w');

sigma = cell2struct(sigma,'Sig');
sigmadt = cell2struct(sigmadt,'Sigdt');

Sig = meanstruct(sigma,'Sig',weights_ub);
Sigdt = meanstruct(sigmadt,'Sigdt',weights_ub);
Sn=meanstruct(sn,'Sn',weights_ub);
Sn.tau=0:(numel(Sn.mean)-1);
