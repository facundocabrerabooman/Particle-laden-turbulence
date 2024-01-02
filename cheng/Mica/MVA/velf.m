function [Velf]=velf(vel,Hess);

% Hess;
length_tresh=512;
j_append=0;

velf=struct();
    velf.good=[];
    velf.lmin=length_tresh;

for j=1:length(vel.data)
    ind=find((vel.data(j).iHes)<Hess);
    if(max(diff(ind))>1)
        indl=find(diff(ind)>1);
        l=diff(indl);
        [m,im]=max(l);
    else
        m=length(ind);
%         indl=[ind(1),ind(length(ind))-1];
        im=1;
    end
    
    if(m > length_tresh)
        velf.good=[vel.good j+j_append];
        ii=ind(indl(im)+1:indl(im+1));
        velf.data(j+j_append).seg=vel.data.seg(ii);
        velf.data(j+j_append).freq=freq(ii);
        velf.data(j+j_append).iHes=iHes(ii);
        velf.i0(j+j_append)=i0(j)+ind(1);
        velf.length(j+j_append)=m;
        velf.file(j+j_append)=cellstr(fname);
        velf.status(j+j_append)=1;
    else
        velf.data(j+j_append).seg=vel.data(j).seg;
        velf.data(j+j_append).freq=vel.data(j).freq;
        velf.data(j+j_append).iHes=vel.data(j).iHes;
        velf.i0(j+j_append)=vel.i0(j);
        %vel.length(j+j_append)=id(j);
%         velf.length(j+j_append)=vel.length(freq);
%         velf.file(j+j_append)=cellstr(fname);
%         velf.status(j+j_append)=0;
    end
end
    