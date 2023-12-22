function velout=vel_trie3(vel)

%%
M=arrayfun(@(x)(mean(x.iHes)),vel.data);
Mf=mean(medfilt1(M,4));
Sf=std(medfilt1(M,4));
HessTh=1;%reject events with mean iHess above 1 sigma
Iout=find(M-Mf>HessTh*Sf);
vel.data(Iout)=[];
vel.i0(Iout)=[];
vel.file(Iout)=[];
%%
M=arrayfun(@(x)(max(abs(x.seg))),vel.data);
ampTh=mean(M)+std(M);%reject events with maximum amplitude below average + 1 sigma
Iout=find(M<ampTh);
vel.data(Iout)=[];
vel.i0(Iout)=[];
vel.file(Iout)=[];
%%
limit_Hess=Mf-.25*Sf;
length_thresh=512;

velout=struct();
n=0;
for j=1:length(vel.data)
    ind=find(vel.data(j).iHes<limit_Hess);
    if(max(diff(ind))>1)
        indl=[find(diff(ind)>1)];
        if(vel.data(j).iHes(ind(numel(ind)))<limit_Hess)
            indl=[indl numel(ind)];
        end
        if (length(indl)==1)
            if(indl<length(ind)/2)
                m=length(ind)-indl;
                indl=[indl,length(ind)];
                im=1;
            else
                m=indl;
                indl=[0,indl];
                im=1;
            end
        else
            indl=[0 indl];
            l=diff(indl);
            [m,im]=max(l);
        end
    else
        m=length(ind);
        indl=[0,length(ind)];
        im=1;
    end
    if(m > length_thresh)
        n=n+1;
        velout.good(n)=n;
        ii=ind(indl(im)+1:indl(im+1));
        velout.data(n).seg=vel.data(j).seg(ii);
        velout.data(n).freq=vel.data(j).freq(ii);
        velout.data(n).iHes=vel.data(j).iHes(ii);
        velout.i0(n)=vel.i0(j)+ii(1);
        velout.length(n)=length(ii);
        velout.file(n)=vel.file(j);
    end
    
    velout.param(3)=ampTh;
    velout.param(2)=limit_Hess;
    
    %     if   ((mod(j,10)==0)||(j==length(vel.data)))
    %          eval(['save ' sname ' vel;']);
    %    end
end

