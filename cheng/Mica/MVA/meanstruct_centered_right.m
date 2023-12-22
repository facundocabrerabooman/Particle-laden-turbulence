function meanstr=meanstruct_centered_right(structarray,field,varargin)

% meanstr=meanstruct_centered(structarray,field,varargin);
%
% calcule, a partir d'un tableau structure structarray(:).field la moyenne
% du champ field centre.
%



if nargin>3
	ind=varargin{2};
	structarray=structarray(ind);
end

ll=arrayfun(@(x)(length(x.(field))),structarray);


Right=arrayfun(@(x)(x.(field)(floor(length(x.(field))/2)+1:numel(x.(field)))),structarray,'UniformOutput',false);
Right=cell2struct(Right',field);
if nargin>2
	fieldweights=varargin{1};
	if ischar(varargin{1})
		Right.(fieldweights)=structarray.(fieldweights);
	end
	Rightm=meanstruct(Right,field,fieldweights);
else
	Rightm=meanstruct(Right,field);
end
%meanstr.N=Leftm.N+Rightm.N;
%meanstr.mean=[Leftm.mean fliplr(Rightm.mean)];
meanstr.N=[Rightm.N];
meanstr.mean=[Rightm.mean] ;
meanstr.std=[Rightm.std] ;