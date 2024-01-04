function pdf=mkpdf5(data,field,n,varargin)

%pdf=mkpdf5(data,field,n,varargin)
%
% computes de PDFs (normalized and not) and statistics from data structure
% where :
% data(i).x is the data
%
%
% pdf.N
% pdf.pdf
% pdf.xpdf
% pdf.pdfn
% pdf.xpdfn
% pdf.mean
% pdf.std
% pdf.skewness
% pdf.flatness


t0=cputime;

%x=zeros(1,data(numel(data)).lc);
%ifin=0;
%for j=1:numel(data)
%        ideb=ifin+1;
%        ifin=ideb+data(j).l-1;
%        x(ideb:ifin)=data(j).x;
        %pour centrer les signaux par la moyenne locale decommenter la
        %ligne ci-dessous
        %x(ideb:ifin)=(x(ideb:ifin)-mean(x(ideb:ifin))); 
     
        %plot(ideb:ifin,x(ideb:ifin));
%end
%moy=mean(x);
%sigma=std(x);

L=arrayfun(@(X)(numel(X.(field))),data);

moy=sum(arrayfun(@(X)(mean(X.(field))),data).*L)/sum(L);
sigma=sqrt(sum(arrayfun(@(X)(mean((X.(field)-moy).^2)),data).*L)/sum(L));



if (nargin>3)&(~isempty(varargin{1}))
    w=varargin{1};
    mm = moy -w*sigma;
    MM = moy +w*sigma;
else
    MM=max(arrayfun(@(X)(max(X.(field))),data));
    mm=min(arrayfun(@(X)(min(X.(field))),data));
    w=round((MM-mm)/sigma);
	%w=round(max(abs(MM),abs(mm))/sigma)
end

%figure;plot(x);
%Xx=linspace(moy-w*sigma,moy+w*sigma,n);
Xx=linspace(mm,MM,n);
%Xx=linspace(moy-w*sigma,moy+w*sigma,n);
Nx=zeros(1,n);

%Lx=0;

 %for jj=1:numel(data)
 %    Nx=Nx+hist(data(jj).(field),Xx);
     %Lx=Lx+numel(data(jj).(field));
 %end
%Lx

NNx=arrayfun(@(X)(hist(X.(field),Xx)),data,'UniformOutput',false);
Nx=sum(cell2mat(NNx'),1);

%size(Nx)
%size(Xx)

pdf.N=L;
pdf.pdf=Nx/integ(Nx,Xx);
pdf.xpdf=Xx;
pdf.meanb=moy;
pdf.mean=integ(pdf.xpdf.*pdf.pdf,pdf.xpdf);
pdf.std=sqrt(integ((pdf.xpdf-pdf.mean).^2.*pdf.pdf,pdf.xpdf));
pdf.stdb=sigma;
pdf.skewness=integ((pdf.xpdf-pdf.mean).^3.*pdf.pdf,pdf.xpdf)/pdf.std^3;
%pdf.skewnessb=skewness(x);
pdf.flatness=integ((pdf.xpdf-pdf.mean).^4.*pdf.pdf,pdf.xpdf)/pdf.std^4;
%pdf.flatnessb=kurtosis(x);

Nx=zeros(1,n);
Xx=linspace(-w,w,n);

%for jj=1:numel(data)
%    Nx=Nx+hist((data(jj).(field)-pdf.mean)/pdf.std,Xx);
%end

NNx=arrayfun(@(X)(hist((X.(field)-moy)/sigma,Xx)),data,'UniformOutput',false);
Nx=sum(cell2mat(NNx'),1);

pdf.pdfn=Nx/integ(Nx,Xx);
pdf.xpdfn=Xx;
%disp(cputime-t0);
   