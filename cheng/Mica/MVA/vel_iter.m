function [vel]=vel_iter(vel,limit_Hess,NbMic,K)

%[velf]=vel_trie(vel,limit_Hess)

length_tresh=512;
vel=vel_trie(vel,limit_Hess);


% dirname='.';
% sname=sprintf('%s/serie%i_vel_N7_K13_%0.5g.mat',dirname,serie,limit_Hess);
good=vel.good;
vel.good=[];
for j=1:length(good)
    [B,A]=butter(4,[min(vel.data(good(j)).freq(K:length(vel.data(good(j)).freq)-K)) max(vel.data(good(j)).freq(K:length(vel.data(good(j)).freq)-K))]*2);
    fdataf=filtfilt(B,A,vel.data(good(j)).seg);
    [freq,iHes]=MVA13mult1(conj(fdataf'),1,NbMic,K);
        
    ind=find(vel.data(good(j)).iHes<limit_Hess);
    if(max(diff(ind))>1)
        indl=[find(diff(ind)>1)];
        if(vel.data(good(j)).iHes(ind(numel(ind)))<limit_Hess)
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
            l=diff(indl);
            [m,im]=max(l);
        end
    else
        m=length(ind);
        indl=[0,length(ind)];
        im=1;
    end
    
    if(m > length_tresh)
        vel.good=[vel.good good(j)];
        ii=ind(indl(im)+1:indl(im+1));
        vel.data(good(j)).seg=vel.data(good(j)).seg(ii);
        vel.data(good(j)).freq=freq(ii);
        vel.data(good(j)).iHes=iHes(ii);
        vel.i0(good(j))=vel.i0(good(j))+ii(1);
        vel.length(good(j))=length(freq(ii));
        vel.file(good(j))=vel.file(good(j));
        vel.status(good(j))=1;
    else
        vel.data(good(j)).seg=vel.data(good(j)).seg;
        vel.data(good(j)).freq=vel.data(good(j)).freq;
        vel.data(good(j)).iHes=vel.data(good(j)).iHes;
        vel.i0(good(j))=vel.i0(good(j));
%         vel.length(j)=vel.id(j);
        vel.length(good(j))=vel.length(good(j));
        vel.file(good(j))=vel.file(good(j));
        vel.status(good(j))=0;
    end
    %     if   ((mod(j,10)==0)||(j==length(vel.data)))
    %          eval(['save ' sname ' vel;']);
    %    end
end

vel=vel_trie(vel,limit_Hess/50);


