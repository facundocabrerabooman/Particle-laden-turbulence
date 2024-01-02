function [M,C]=cumulants_2(inc,p);

% [M]=cumulants_2(inc,p);
%
% calcule les moments d'ordre p des incréments à partir de leur PDF

M=zeros(numel(p),numel(inc.dt));

for kk=1:numel(inc.dt)
    NN=integ(inc.pdf(kk).pdf,inc.pdf(kk).xpdf);
    for ii=1:numel(p)
        M(ii,kk)=2*integ(inc.pdf(kk).pdfn.*inc.pdf(kk).xpdfn.^p(ii),inc.pdf(kk).xpdfn); 
    end

    P=polyfit(log(M(:,kk))',p,3);
    C.c0(kk)=P(4);
    C.c1(kk)=P(3);
    C.c2(kk)=P(2);
    C.c3(kk)=P(1);
end


