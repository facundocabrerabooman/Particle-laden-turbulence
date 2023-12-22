function meanstr=meanstruct_max(structarray,field,varargin)

% meanstr=meanstruct_centered(structarray,field,varargin);
%
% calcule, a partir d'un tableau structure structarray(:).field la moyenne
% du champ field centre sur le max.
%



if nargin>3
	ind=varargin{2};
	structarray=structarray(ind);
end



ll=arrayfun(@(x)(length(x.(field))),structarray);
[m,I]=arrayfun(@(x)(max(x.(field))),structarray)

Left=arrayfun(@(x,J)(flipud(fliplr(x.(field)(1:J)))),structarray,I,'UniformOutput',false);
Left=cell2struct(Left,field);

if nargin>2
	fieldweights=varargin{1};
	Left.(fieldweights)=structarray.(fieldweights);
	Leftm=meanstruct(Left,field,varargin);
else
	Leftm=meanstruct(Left,field);
end

Right=arrayfun(@(x,J)(x.(field)(J+1:numel(x.(field)))),structarray,I,'UniformOutput',false);
Right=cell2struct(Right,field);
if nargin>2
	fieldweights=varargin{1};
	Right.(fieldweights)=structarray.(fieldweights);
	Rightm=meanstruct(Right,field,varargin);
else
	Rightm=meanstruct(Right,field);
end
%meanstr.N=Leftm.N+Rightm.N;
%meanstr.mean=[Leftm.mean fliplr(Rightm.mean)];
meanstr.N=[fliplr(flipud(Leftm.N)) Rightm.N];
meanstr.mean=[fliplr(flipud(Leftm.mean)) Rightm.mean] ;
