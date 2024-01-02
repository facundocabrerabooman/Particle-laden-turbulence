function plotLagr(LagFilename,mycolormap,Fs,n,if2layers,ifsavefig)

load(LagFilename)

fLagr = './Figures_Lag';
if ifsavefig==1
    mkdir(fLagr)
end

%%
plotLagr_std(w,s,wopt,mycolormap,ifsavefig,fLagr)

%% normalized pdfs
plotLagr_pdfs(LagragianStats.pdfV,LagragianStats.pdfA,mycolormap,ifsavefig,fLagr)
% plotLagr_pdfs_lavision(pdfV,pdfA,mycolormap1,1) % compare to Lavision output

%%
plotLagr_MSDs(LagragianStats.MSD,Fs,mycolormap,ifsavefig,fLagr)

%%
plotLagr_S2L(LagragianStats.S2L,Fs,mycolormap,ifsavefig,fLagr)

%% Correlation function -- using Gaussian filter
plotLagr_Corr(LagragianStats.Ruu,LagragianStats.Raa,Fs,n,if2layers,mycolormap,ifsavefig,fLagr)

%% Correlation function -- using dt method (denosied)
plotLagr_dtCorr(LagragianStats.tau,LagragianStats.corrv,LagragianStats.corra,ifsavefig,fLagr)

