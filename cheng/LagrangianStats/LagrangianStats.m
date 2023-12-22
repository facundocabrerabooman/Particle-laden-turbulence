function LagragianStats = LagrangianStats(track,Fs,Nbins,n,ndts,ndtl,ifsave)

%% 1 time - 1 particle statistics
% velocity and acceleration pdfs
disp('calculating PDFs')
[LagragianStats.pdfV,LagragianStats.pdfA] = LagrangianStatPDF(track,Nbins,n);

%% 2 times - 1 particle statistics (Lagrangian statistics)
%% MSD
disp('calculating MSDs')
LagragianStats.MSD = LagrangianStatMSD(track);

%% 2nd SF -- Lagrangian
disp('calculating S2L')
LagragianStats.S2L = LagrangianStatS2L(track);

%% Correlation
disp('calculating Correlation functions')
[LagragianStats.Ruu,LagragianStats.Raa] = LagrangianStatCorr(track);

%%  Correlation function -- using dt method (denosied, Romain)
disp('calculating Correlation functions -- dt method')
LagragianStats.ndts = ndts;
LagragianStats.ndtl = ndtl;
[LagragianStats.tau,LagragianStats.corrv,LagragianStats.corra] = dtCorr(track,ndts,ndtl,Fs);

%% Lagrangian data&plots
LagFilename = './LagrangianStats.mat';
if ifsave==1
    disp('Saving')
    save(LagFilename,'LagragianStats','-append')
end

