function [Sn Sig Sigdt]=structFunc_struct(S1,field1,n,varargin)
%
% Sn=SructFunc(S1,field1,n)
%

if nargin>3
    Nmax = varargin{1};
else
    N = arrayfun(@(X)(numel(X.(field1))),S1);
    Nmax = numel(N);
end


disp(sprintf('Moyennes bulle a bulle ...\n'));
[sn, weights_ub sigma sigmadt]=arrayfun(@(x)(structFunc(x.(field1)(1:min(Nmax,end)),n)),S1,'UniformOutput',false);

disp(sprintf('Moyenne d''ensemble ...\n'));
sn=cell2struct(sn,'Sn');
weights_ub=cell2struct(weights_ub,'w');

sigma = cell2struct(sigma,'Sig');
sigmadt = cell2struct(sigmadt,'Sigdt');

Sig = meanstruct(sigma,'Sig',weights_ub);
Sigdt = meanstruct(sigmadt,'Sigdt',weights_ub);
Sn=meanstruct(sn,'Sn',weights_ub);
Sn.tau=0:(numel(Sn.mean)-1);

