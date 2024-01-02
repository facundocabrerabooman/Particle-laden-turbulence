function meanstr=meanstruct(structarray,field,varargin)

% meanstruct=meanstruct(structarray,field,varargin);
%
% calcule, a partir d'un tableau structure structarray(:).field la moyenne
% du champ field.
%

flagweights=0;
if nargin>2
	if ischar(varargin{1})
		fieldweights=varargin{1};
		flagweights=1;
	else
		weights=varargin{1};
		if ~isempty(weights)
			flagweights=1;
		end
	end
end

ll=arrayfun(@(x)(length(x.(field))),structarray);
%[ll,II]=sort(ll);
%structarray=structarray(II);

if nargin>3
	MM=varargin{2};
	if MM==0
		MM=min(ll);
	end
else 
	MM=max(ll);
end

mm=min(ll);

if flagweights==1
	if ~ischar(varargin{1})
		Npts=sum(cell2mat(arrayfun(@(y)(y.w(1:mm)'),weights,'UniformOutput',false)),1);
		meanstr.mean(1:mm)=sum(cell2mat(arrayfun(@(x,y)(x.(field)(1:mm).*y.w(1:mm)'),structarray,weights,'UniformOutput',false)))./Npts;
		meanstr.std(1:mm)=sqrt(sum(cell2mat(arrayfun(@(x,y)((x.(field)(1:mm)-meanstr.mean(1:mm)).^2.*y.w(1:mm)'),structarray,weights,'UniformOutput',false)))./Npts);
		meanstr.l(1:mm)=Npts;
	else
		Npts=sum(cell2mat(arrayfun(@(x)(x.(fieldweights)(1:mm)),structarray,'UniformOutput',false)));
		meanstr.mean(1:mm)=sum(cell2mat(arrayfun(@(x)(x.(field)(1:mm).*x.(fieldweights)(1:mm)'),structarray,'UniformOutput',false)))./Npts;
		meanstr.std(1:mm)=sqrt(cell2mat(sum(arrayfun(@(x)((x.(field)(1:mm)-meanstr.mean(1:mm)).^2.*x.(fieldweights)(1:mm)'),structarray,'UniformOutput',false)))./Npts);
		meanstr.l(1:mm)=Npts;
	end
else
	meanstr.mean(1:mm)=sum(cell2mat(arrayfun(@(x)(x.(field)(1:mm)),structarray,'UniformOutput',false)))/numel(structarray);
	meanstr.std(1:mm)=sqrt(sum(cell2mat(arrayfun(@(x)((x.(field)(1:mm)-meanstr.mean(1:mm)).^2),structarray,'UniformOutput',false)))/numel(structarray));
end
	meanstr.N(1:mm)=numel(structarray);

for kk=mm+1:MM
	J=find(ll-kk>=0);
	meanstr.N(kk)=numel(J);
	if flagweights==1
		if ~ischar(varargin{1})
			Npts=sum(arrayfun(@(y)(y.w(kk)),weights(J)));
			meanstr.mean(kk)=sum(arrayfun(@(x,y)(x.(field)(kk)*y.w(kk)),structarray(J),weights(J)))/Npts;
			meanstr.std(kk)=sqrt(sum(arrayfun(@(x,y)((x.(field)(kk)-meanstr.mean(kk))^2*y.w(kk)),structarray(J),weights(J)))/Npts);
			meanstr.l(kk)=Npts;
		else
			Npts=sum(arrayfun(@(x)(x.(fieldweights)(kk)),structarray(J)));
			meanstr.mean(kk)=sum(arrayfun(@(x)(x.(field)(kk)*x.(fieldweights)(kk)),structarray(J)))/Npts;
			meanstr.std(kk)=sqrt(sum(arrayfun(@(x)((x.(field)(kk)-meanstr.mean(kk))^2*x.(fieldweights)(kk)),structarray(J)))/Npts);
			meanstr.l(kk)=Npts;
		end
	else
		meanstr.mean(kk)=sum(arrayfun(@(x)(x.(field)(kk)),structarray(J)))/numel(J);
		meanstr.std(kk)=sqrt(sum(arrayfun(@(x)((x.(field)(kk)-meanstr.mean(kk))^2),structarray(J)))/numel(J));
	end
end
