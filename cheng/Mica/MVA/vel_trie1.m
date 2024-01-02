function [vel]=vel_trie(vel,limit_Hess)
%[velf]=vel_trie(vel,limit_Hess)

length_tresh=512;
vel.good=[];

% dirname='.';
% sname=sprintf('%s/serie%i_vel_N7_K13_%0.5g.mat',dirname,serie,limit_Hess);

for j=1:length(vel.data)
    ind=find(vel.data(j).iHes<limit_Hess);
    if(max(diff(ind))>1)
        indl=find(diff(ind)>1);
        l=diff(indl);
        [m,im]=max(l);
    else
        m=length(ind);
        indl=[ind(1),ind(length(ind))-1];
        im=1;
    end

    if(m > length_tresh)
        vel.good=[vel.good j];
%         ii=ind(indl(im)+1:indl(im+1));
%         vel.data(j).seg=vel.data(j).seg(ii);
%         vel.data(j).freq=vel.data(j).freq(ii);
%         vel.data(j).iHes=vel.data(j).iHes(ii);
%         vel.i0(j)=vel.i0(j)+ind(1);
%         vel.length(j)=m;
%         vel.file(j)=cellstr(fname);
%         vel.status(j)=1;
%     else
%         vel.data(j).seg=vel.data(j).seg;
%         vel.data(j).freq=vel.data(j).freq;
%         vel.data(j).iHes=vel.data(j).iHes;
%         vel.i0(j)=vel.i0(j);
% %         vel.length(j)=vel.id(j);
%         vel.length(j)=vel.length(j);
%         vel.file(j)=vel.file(j);
%         vel.status(j)=0;
    end
    %     if   ((mod(j,10)==0)||(j==length(vel.data)))
    %          eval(['save ' sname ' vel;']);
    %    end
end

