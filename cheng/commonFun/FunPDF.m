function [ipdf,bincenter] = FunPDF(series,nbins)

[ihist,bincenter] = hist(series,nbins);
ihist_sum = sum(ihist);
dbin = bincenter(2)-bincenter(1);
ipdf = ihist/ihist_sum/dbin;

ipdf = ipdf';
bincenter = bincenter';