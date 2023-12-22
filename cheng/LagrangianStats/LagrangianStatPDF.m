function [pdfV,pdfA] = LagrangianStatPDF(track,Nbins,w)

%% Calculate & plot velocity and acceleration pdfs
pdfV(1) = mkpdf5(track,'Vx',Nbins,w);
pdfV(2) = mkpdf5(track,'Vy',Nbins,w);
pdfV(3) = mkpdf5(track,'Vz',Nbins,w);

pdfA(1) = mkpdf5(track,'Ax',Nbins,2*w);
pdfA(2) = mkpdf5(track,'Ay',Nbins,2*w);
pdfA(3) = mkpdf5(track,'Az',Nbins,2*w);