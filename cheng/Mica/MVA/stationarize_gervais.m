function Sout=stationarize_gervais(S1,field1,moysup,sigsup);

% Sout=stationarize_gervais(S1,field1,Icenter,Profil);
%
% stationarise la vitesse a la P.Gervais : (v-v>)/sig>

ll=arrayfun(@(x)(numel(x.(field1))),S1);

Sout=arrayfun(@(x,L)((x.(field1)-moysup(L))/sigsup(L)),S1,ll,'UniformOutput',false);



