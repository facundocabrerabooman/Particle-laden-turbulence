function pdf=mkpdf4c(data,n,varargin)

%pdf=mkpdf4b(data,n,varargin)
%
% computes de PDFs (normalized and not) and statistics from data structure
% where :
%
%
% pdf.N
% pdf.pdf
% pdf.xpdf
% pdf.mean
% pdf.std
% pdf.skewness
% pdf.flatness


t0=cputime;


L=numel(data);

moy=nanmean(data);
sigma=nanstd(data);

if nargin>2
    type = varargin{1};
    
    switch type
        case ''
        case 'c'
            data = data-moy;
        case 'cr'
            data = (data-moy)/sigma;
        case 'r'
            data = data/sigma;
        case 'n'
            data = data/moy;
    end
end

if nargin>3
    a=varargin{1};
    b=varargin{2};
    Xx=linspace(a,b,n);
else
    Xx=linspace(min(data),max(data),n);
    %Xx=logspace(log10(min(data)+1e-30),log10(max(data)),n);
end

Nx=zeros(1,n);

Nx=hist(data,Xx);
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

