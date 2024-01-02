function meanstr=meanstruct_centered(structarray,field,varargin)

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
Left=arrayfun(@(x)(flipud(fliplr(x.(field)(1:floor(length(x.(field))/2))))),structarray,'UniformOutput',false);
Left=cell2struct(Left',field);
if nargin>2
	fieldweights=varargin{1};
	if ischar(varargin{1})
		Left.(fieldweights)=structarray.(fieldweights);
	end
	Leftm=meanstruct(Left,field,fieldweights);
else
	Leftm=meanstruct(Left,field);
end

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
meanstr.Center=numel(Leftm)+1;
meanstr.N=[flipud(fliplr(Leftm.N)) Rightm.N];
meanstr.mean=[flipud(fliplr(Leftm.mean)) Rightm.mean] ;
meanstr.std=[flipud(fliplr(Leftm.std)) Rightm.std] ;