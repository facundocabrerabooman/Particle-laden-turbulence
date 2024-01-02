function [sigma_a, sigma_v, mean_a, mean_v, w]=findFilterWidth(vel);

% [sigma_a, sigma_v, mean_a, mean_v, w]=findFilterWidth(vel);

w=1:5:30;
l=max(25,5*w);

L=arrayfun(@(x)(numel(x.vit)),vel.velocity);

for j=1:numel(w)
    disp(sprintf('w = %i',w(j)));
    kerp = posfiltcoef(w(j),l(j));
    kerv = velfiltcoef(w(j),l(j));
    vitf=[];
    acc=[];
    for jj=1:numel(L)
        if L(jj)>l(j)+1
            vitftmp=conv(kerp,vel.velocity(jj).vit);
            acctmp=conv(kerv,vel.velocity(jj).vit);
        
            vitf=[vitf vitftmp(l(j):length(vitftmp)-l(j))];
            acc=[acc acctmp(l(j):length(acctmp)-l(j))];
            clear vitftmp acctmp;
        end
    end
    
    sigma_v(j)=std(vitf);
    mean_v(j)=mean(vitf);
    sigma_a(j)=std(acc);
    mean_a(j)=mean(acc);
end

