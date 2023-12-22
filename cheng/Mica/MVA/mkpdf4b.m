function pdf=mkpdf4b(data,n,varargin)

%pdf=mkpdf4b(data,n,varargin)
%
% computes de PDFs (normalized and not) and statistics from data structure
% where :
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

L=numel(data);

moy=mean(data);
sigma=std(data);



if nargin>3
    a=varargin{1};
    b=varargin{2};
    Xx=linspace(a,b,n);
elseif nargin>2
    w=varargin{1};
    Xx=linspace(moy-w*sigma,moy+w*sigma,n);
else
    MM=max(data-moy);
    mm=min(data-moy);
    w=round((MM-mm)/sigma);
    Xx=linspace(min(data),max(data),n);
    %Xx=logspace(log10(min(data)+1e-30),log10(max(data)),n);
end

%figure;plot(x);

Nx=zeros(1,n);
%Lx=0;

 %for jj=1:numel(data)
 %    Nx=Nx+hist(data(jj).x,Xx);
     %Lx=Lx+numel(data(jj).x);
 %end
%Lx

NNx=hist(data,Xx);
%Nx=sum(NNx,1);
Nx=NNx;
dXx=diff(Xx);
dXx(end+1)=dXx(end);

pdf.N=L;
%pdf.pdf=Nx/integ(Nx,Xx);
pdf.pdf=Nx./dXx/sum(Nx);
pdf.xpdf=Xx;
%pdf.meanb=moy;
pdf.mean=integ(pdf.xpdf.*pdf.pdf,pdf.xpdf);
pdf.std=sqrt(integ((pdf.xpdf-pdf.mean).^2.*pdf.pdf,pdf.xpdf));
%pdf.stdb=sigma;
pdf.skewness=integ((pdf.xpdf-pdf.mean).^3.*pdf.pdf,pdf.xpdf)/pdf.std^3;
%pdf.skewnessb=skewness(x);
pdf.flatness=integ((pdf.xpdf-pdf.mean).^4.*pdf.pdf,pdf.xpdf)/pdf.std^4;
%pdf.flatnessb=kurtosis(x);

Nx=zeros(1,n);
if nargin>3
Xx=linspace(a,b,n);    
else
Xx=linspace(-w,w,n);
end

%for jj=1:numel(data)
%    Nx=Nx+hist((data(jj).x-pdf.mean)/pdf.std,Xx);
%end
NNx=hist((data-moy)/sigma,Xx);
%Nx=sum(cell2mat(NNx),1);
Nx=NNx;

pdf.pdfn=Nx/integ(Nx,Xx);
pdf.xpdfn=Xx;
%disp(cputime-t0);
   