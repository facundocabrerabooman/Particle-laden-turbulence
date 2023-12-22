function [i,vali,delta]=findi(var,val,quant)

%
% [i,vali,delta]=findi(var,val,quant)
%
% cette fonction est un "inverseur" indice-valeur
%
% elle retourne une liste de taille 'quant' des indices i tels que var(i)~val
% vali renvoie les valeurs var(i)
% delta renvoie les valeurs de |var(i)-val|
%

[Y,ind]=sort(var);
i=[];
j=1;
n=0;
while(n<quant);
    n=n+1;
    i1=j-1+min(find(Y(j:length(Y))>=val));
    i2=j-1+max(find(Y(j:length(Y))<=val));
    if (~isempty(i1)&~isempty(i2))
        if abs(Y(i1)-val)<abs(Y(i2)-val)
            i=[i1,i];
        else
            i=[i2,i];
        end
    else
        if isempty(i1),i=[i2,i];end;
        if isempty(i2),i=[i1,i];end;
    end
    j=i(n)+1;
end
i=ind(i);
vali=var(i);
delta=abs(var(i)-val);


    
