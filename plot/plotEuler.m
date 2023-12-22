function plotEuler(eulerStats,Fs,mycolormap,ifsubsmean,ifsavefig)


if ifsubsmean == 1
	fEuler = './Figures_Euler_subsMean';
elseif ifsubsmean == 0
	fEuler = './Figures_Euler';
end
mkdir(fEuler)
%% check stationary
plotEuler_stationary(eulerStats,Fs,mycolormap,ifsavefig,fEuler)

%% S2x, S2y, S2z
plotEuler_S2E(eulerStats,mycolormap,ifsavefig,fEuler)

%% Ruu
plotEuler_Ruur(eulerStats,mycolormap,ifsavefig,fEuler)

%% PSD
plotEuler_PSDs(eulerStats,mycolormap,ifsavefig,fEuler)

%% Sau & Saulong
plotEuler_Sau(eulerStats,mycolormap,ifsavefig,fEuler)

%% Splong
plotEuler_Spn(eulerStats,mycolormap,ifsavefig,fEuler)

%% SplongAbs
plotEuler_SpnAbs(eulerStats,mycolormap,ifsavefig,fEuler)

%% compensate Eulerian plots
plotEuler_epsilon(eulerStats,mycolormap,ifsavefig,fEuler)