function [vel]=vel_trie(vel,limit_Hess,amp_thresh)
%[velf]=vel_trie(vel,limit_Hess)

length_thresh=512;
%amp_thresh=1e-8;
vel.good=[];

% dirname='.';
% sname=sprintf('%s/serie%i_vel_N7_K13_%0.5g.mat',dirname,serie,limit_Hess);

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
    if((m > length_thresh)&(max(abs(vel.data(j).seg).^2)>amp_thresh))
        vel.good=[vel.good j];
        ii=ind(indl(im)+1:indl(im+1));
        vel.data(j).seg=vel.data(j).seg(ii);
        vel.data(j).freq=vel.data(j).freq(ii);
        vel.data(j).iHes=vel.data(j).iHes(ii);
        vel.i0(j)=vel.i0(j)+ii(1);
        vel.length(j)=length(ii);
        vel.file(j)=vel.file(j);
        vel.status(j)=1;
    else
        vel.data(j).seg=vel.data(j).seg;
        vel.data(j).freq=vel.data(j).freq;
        vel.data(j).iHes=vel.data(j).iHes;
        vel.i0(j)=vel.i0(j);
        %vel.length(j)=vel.id(j);
        vel.length(j)=vel.length(j);
        vel.file(j)=vel.file(j);
        vel.status(j)=0;
    end
    
    vel.param(3)=amp_thresh;
    vel.param(2)=limit_Hess;
    
    %     if   ((mod(j,10)==0)||(j==length(vel.data)))
    %          eval(['save ' sname ' vel;']);
    %    end
end

