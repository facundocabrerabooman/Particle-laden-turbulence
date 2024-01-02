function meanstr=meanstruct_centered2(structarray,field,Icenter,varargin)

% meanstr=meanstruct_centered(structarray,field,varargin);
%
% calcule, a partir d'un tableau structure structarray(:).field la moyenne
% du champ field centre autour de Icenter;
%



if nargin>4
	ind=varargin{2};
	structarray=structarray(ind);
end

ll=arrayfun(@(x)(length(x.(field))),structarray);
Left=arrayfun(@(x,I)(flipud(fliplr(x.(field)(1:I)))),structarray,Icenter,'UniformOutput',false);
Left=cell2struct(Left,field);

max(arrayfun(@(x)(numel(x.(field))),Left))


if nargin>3
	fieldweights=varargin{1};
	Left.(fieldweights)=structarray.(fieldweights);
	Leftm=meanstruct(Left,field,varargin);
else
	Leftm=meanstruct(Left,field);
end

Right=arrayfun(@(x,I)(x.(field)(I+1:numel(x.(field)))),structarray,Icenter,'UniformOutput',false);
Right=cell2struct(Right,field);
if nargin>3
	fieldweights=varargin{1};
	Right.(fieldweights)=structarray.(fieldweights);
	Rightm=meanstruct(Right,field,varargin);
else
	Rightm=meanstruct(Right,field);
end
%meanstr.N=Leftm.N+Rightm.N;
%meanstr.mean=[Leftm.mean fliplr(Rightm.mean)];
meanstr.Center=numel(Leftm)+1;
meanstr.N=[flipud(fliplr(Leftm.N)) Rightm.N];
meanstr.mean=[flipud(fliplr(Leftm.mean)) Rightm.mean] ;
meanstr.std=[flipud(fliplr(Leftm.std)) Rightm.std] ;
