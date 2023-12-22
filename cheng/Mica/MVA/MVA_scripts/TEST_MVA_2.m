% load DATA
load Sig_Bulle_4

% MVA Mordant Pinton Michel
sig = conj(diff(detrend(Sig_Bulle_4_i_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.1;
wMVA;
Dop_MVA_Bulle_4_i = freq(1,:)*F_Ech;
InvHess_MVA_Bulle_4_i = iHes;

sig = conj(diff(detrend(Sig_Bulle_4_j_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.05;
wMVA;
Dop_MVA_Bulle_4_j = freq(1,:)*F_Ech;
InvHess_MVA_Bulle_4_j = iHes;

sig = conj(diff(detrend(Sig_Bulle_4_k_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.05;
wMVA;
Dop_MVA_Bulle_4_k = freq(1,:)*F_Ech;
InvHess_MVA_Bulle_4_k = iHes;

sig = conj(diff(detrend(Sig_Bulle_4_l_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.1;
wMVA;
Dop_MVA_Bulle_4_l = freq(1,:)*F_Ech;
InvHess_MVA_Bulle_4_l = iHes;
clear sig

REV_MVA_FLG = 0;
if (REV_MVA_FLG)
	sig = flipud(diff(detrend(Sig_Bulle_4_i_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_i = fliplr(freq(1,:))*F_Ech;
	InvHess_REV_MVA_Bulle_4_i = fliplr(iHes);
	
	sig = flipud(diff(detrend(Sig_Bulle_4_j_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_j = fliplr(freq(1,:))*F_Ech;
	InvHess_REV_MVA_Bulle_4_j = fliplr(iHes);
	
	sig = flipud(diff(detrend(Sig_Bulle_4_k_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_k = fliplr(freq(1,:))*F_Ech;
	InvHess_REV_MVA_Bulle_4_k = fliplr(iHes);
	
	sig = flipud(diff(detrend(Sig_Bulle_4_l_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_l = fliplr(freq(1,:))*F_Ech;
	InvHess_REV_MVA_Bulle_4_l = fliplr(iHes);
end
figure
MAX_DOP_FREQ = max([Dop_MVA_Bulle_4_i Dop_MVA_Bulle_4_j Dop_MVA_Bulle_4_k Dop_MVA_Bulle_4_l]);
MAX_DOP_FREQ = ceil(MAX_DOP_FREQ/1000)*1000;
subplot(221); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_MVA_Bulle_4_i,'b'); axis([0 0.125 0 MAX_DOP_FREQ]); grid
subplot(222); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_MVA_Bulle_4_j,'g'); axis([0 0.125 0 MAX_DOP_FREQ]); grid
subplot(223); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_MVA_Bulle_4_k,'r'); axis([0 0.125 0 MAX_DOP_FREQ]); grid
subplot(224); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_MVA_Bulle_4_l,'k'); axis([0 0.125 0 MAX_DOP_FREQ]); grid

% Clean MVA variables
clear P MaxIter NbMic NbS Nb_Iter PiB PiBc Puis Reps Rx Ry Sc Teta Tetac U V VP Veps epsilon freq iHes iHess iSSt ii k l ll progress
clear sig syg syg2 tot vp window Gamma Gammac Grad Hess0 Hessian A DSc K MVAinit

% Extraction des 3 composantes xyz de vitesse
Psi = 26*pi/180; Phi = atan(sin(Psi)); Cson = 340;
F_inc_i = 110000; F_inc_j = 122000; F_inc_k = 110000; F_inc_l = 122000;

Dop2Vit_i = Cson/(2*sin(Phi))/F_inc_i;
Dop2Vit_j = Cson/(2*sin(Phi))/F_inc_j;
Dop2Vit_k = Cson/(2*sin(Phi))/F_inc_k;
Dop2Vit_l = Cson/(2*sin(Phi))/F_inc_l;

Vx =   (Dop2Vit_l*Dop_MVA_Bulle_4_l - Dop2Vit_k*Dop_MVA_Bulle_4_k)/(2*sin(Psi));
Vy =   (Dop2Vit_i*Dop_MVA_Bulle_4_i - Dop2Vit_j*Dop_MVA_Bulle_4_j)/(2*sin(Psi));
Vz_1 = (Dop2Vit_l*Dop_MVA_Bulle_4_l + Dop2Vit_k*Dop_MVA_Bulle_4_k)/(2*cos(Psi));
Vz_2 = (Dop2Vit_i*Dop_MVA_Bulle_4_i + Dop2Vit_j*Dop_MVA_Bulle_4_j)/(2*cos(Psi));

% Overlay REV_MVA
if (REV_MVA_FLG),
	subplot(221); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_REV_MVA_Bulle_4_i,'b-.'); hold off
	subplot(222); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_REV_MVA_Bulle_4_j,'g-.'); hold off
	subplot(223); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_REV_MVA_Bulle_4_k,'r-.'); hold off
	subplot(224); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech,Dop_REV_MVA_Bulle_4_l,'k-.'); hold off
end
