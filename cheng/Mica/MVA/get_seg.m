function [segments,i0,id,im]=get_seg(s,varargin)

%segments=get_seg(s);

min_length=2000;
deborde_droite=0;	%etend la longueur des segments trouver d'autant de part et d'autre
deborde_gauche=0;
Ngroupe=1;		%regroupe Ngroupe segments en 1;

if nargin>1
	s2=varargin{1};
else
	s2=s;
end


%s2=s;
%s2=abs(s).^2;
%figure;plot(s2);

thresh1=median(s2);
thresh2=2.*max(s2);



%medfDdata=1000.*median(s2);
medfDdata=thresh2/10
%medfDdata=0;

%[B,A]=butter(4,.005);
%s2=filtfilt(B,A,s2);

ns=0;
    
I=(2:length(s)-1);
sp=s2;
sp(I)=sp(I+1);
id0m=find(s2<medfDdata&sp>medfDdata); %(points passant par zero ou s est croissant)
id0p=find(s2>medfDdata&sp<medfDdata); %(points passant par zero ou s est decroissant)

if numel(id0p)<numel(id0m)
    id0p(numel(id0p)+1)=numel(s2);
end

if Ngroupe~=1
	id0m=id0m(1:Ngroupe:numel(id0m));
	id0p=id0p(1:Ngroupe:numel(id0p));
end

if (numel(id0p)*numel(id0m)==0)
    segments=[];
    i0=[];
    id=[];
    xc=[];
    return;
end

if id0p(1)<id0m(1)
    I=(1:length(id0p)-1);
    id0p(I)=id0p(I+1);
end

if s(length(s))<0;
    iend=length(id0m)-1;
else
    iend=length(id0m);
end

for i=1:iend
    
    [maxs imaxs]=max(s2(id0m(i):id0p(i)));
	
    if  maxs>=thresh1&maxs<=thresh2
        idd=id0p(i)-id0m(i)+1;
        %idxmaxs=[idxmaxs id];
        
    	if(idd>min_length)
            
           itop=find(s2(id0m(i):id0p(i))>maxs/5);
            ns=ns+1;
            %xc(ns)=xcontrast(abs(s(id0m(i)+itop-1)).^2./s2(id0m(i)+itop-1));
			i0(ns)=max(1,id0m(i)-deborde_gauche);
			i0p=min(numel(s),id0p(i)+deborde_droite);
            id(ns)=i0p-i0(ns)+1;
			im(ns)=id0m(i)+imaxs-1;
			segments(ns).seg=s(i0(ns):i0p);
            
            
            %segments(ns).i0=id0m(i);
            %segments(ns).length=id;
        end                    
    end
end

if ~exist('segments','var')
    segments=[];
    i0=[];
    id=[];
    xc=[];
end
    


function xc=xcontrast(f1)
xc=abs((max(f1)-min(f1))/(max(f1)+min(f1)));
