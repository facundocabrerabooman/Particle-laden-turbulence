function Sout=minus_centered(S1,field1,Icenter,B,BCentre);

ll=arrayfun(@(x)(numel(x.(field1))),S1);

Ideb=BCentre-Icenter;
Ifin=BCentre+ll-Icenter-1;



Sout=arrayfun(@(x,Id,If)(x.(field1)-B(Id:If)),S1,Ideb,Ifin,'UniformOutput',false);

Sout=cell2struct(Sout,field1);

