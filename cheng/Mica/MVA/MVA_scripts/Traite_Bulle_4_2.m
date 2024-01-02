% Bulle 4, Decimate x4
% Bulle_4 = 103000; Pos = Bulle_4; s_Bulle_4_dec4 = Pos-2048:Pos+2048;
% Sig_Bulle_4_i_dec4 = Voie_i_dec4(s_Bulle_4_dec4);
% Sig_Bulle_4_j_dec4 = Voie_j_dec4(s_Bulle_4_dec4);
% Sig_Bulle_4_k_dec4 = Voie_k_dec4(s_Bulle_4_dec4);
% Sig_Bulle_4_l_dec4 = Voie_l_dec4(s_Bulle_4_dec4);
load Sig_Bulle_4

% Test Temps-Frequence
% Amplitudes
figure
subplot(421);plot((0:length(s_Bulle_4_dec4)-2)/F_Ech*4, abs(diff(Sig_Bulle_4_i_dec4)).^2); grid
title(sprintf('Bulle #4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(423);plot((0:length(s_Bulle_4_dec4)-2)/F_Ech*4, abs(diff(Sig_Bulle_4_j_dec4)).^2); grid
subplot(425);plot((0:length(s_Bulle_4_dec4)-2)/F_Ech*4, abs(diff(Sig_Bulle_4_k_dec4)).^2); grid
subplot(427);plot((0:length(s_Bulle_4_dec4)-2)/F_Ech*4, abs(diff(Sig_Bulle_4_l_dec4)).^2); grid
% Spectres
F_tmp = (-(length(s_Bulle_4_dec4)-1)/2:(length(s_Bulle_4_dec4)-1)/2-1)/(length(s_Bulle_4_dec4)-1);
Spect_Bulle_4_i_dec4 = (abs(fftshift(fft(diff(detrend(conj(Sig_Bulle_4_i_dec4))).*hanning(length(s_Bulle_4_dec4)-1))))).^2;
Spect_Bulle_4_j_dec4 = (abs(fftshift(fft(diff(detrend(conj(Sig_Bulle_4_j_dec4))).*hanning(length(s_Bulle_4_dec4)-1))))).^2;
Spect_Bulle_4_k_dec4 = (abs(fftshift(fft(diff(detrend(conj(Sig_Bulle_4_k_dec4))).*hanning(length(s_Bulle_4_dec4)-1))))).^2;
Spect_Bulle_4_l_dec4 = (abs(fftshift(fft(diff(detrend(conj(Sig_Bulle_4_l_dec4))).*hanning(length(s_Bulle_4_dec4)-1))))).^2;
subplot(422);plot(F_tmp*F_Ech/4, 10*log10(Spect_Bulle_4_i_dec4)); axis([-5000 5000 -60 -10]); grid
title(sprintf('Bulle #4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(424);plot(F_tmp*F_Ech/4, 10*log10(Spect_Bulle_4_j_dec4)); axis([-5000 5000 -60 -10]); grid
subplot(426);plot(F_tmp*F_Ech/4, 10*log10(Spect_Bulle_4_k_dec4)); axis([-5000 5000 -60 -10]); grid
subplot(428);plot(F_tmp*F_Ech/4, 10*log10(Spect_Bulle_4_l_dec4)); axis([-5000 5000 -60 -10]); grid

% CLEAN
CLEAN_FLG = 1
if (CLEAN_FLG)
	clear Voie_i Voie_j Voie_k Voie_l data
end

% TFR
% Test tfr
% Time_test = (0:length(Sig_Bulle_4_i_dec4)-1);
% Sig_Test = fmlin(length(Time_test), 0.1, 0.4);
% h = hanning(127);
% [tfr_test, rtfr_test]=  mex_tfrrsp(Sig_Test, 1:length(Sig_Bulle_4_i_dec4), 512, h, 1);
% h=hanning(15);
% figure
% plot(F_Bulles/F_Ech*4, 10*log10(P_bulles_i),'b')
% F_Dop_Avg = 0.08; F_Dop_test = F_Dop_Avg + 0.05*sin(2*pi*Time_test*0.25);
% Sig_test = exp(i*2*pi*Time_test.*F_Dop_test);
% Spect_Sig_Test = (abs(fftshift(fft(Sig_test.*hanning(length(s_Bulle_4_dec4)-1))))).^2;
% plot(
% sig_test = 

% [tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_i_dec4)), 1:4096, 1024, hanning(513), 0);
% RTFR1024_Bulle_4_han513_i = rtfr_tmp(1:250,:);
% TFR1024_Bulle_4_han513_i = tfr_tmp(1:250,:);
% [tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_i_dec4)), 1:4096, 1024, hanning(255), 0);
% RTFR1024_Bulle_4_han255_i = rtfr_tmp(1:250,:);
% TFR1024_Bulle_4_han255_i = tfr_tmp(1:250,:);

% TFRSPet RTFRSP sur 1024 points

%[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_i_dec4)), 1:4096, 1024, hanning(127), 0);
%TFR1024_Bulle_4_han127_i = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_i = rtfr_tmp(1:250,:); 
%[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_j_dec4)), 1:4096, 1024, hanning(127), 0);
%TFR1024_Bulle_4_han127_j = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_j = rtfr_tmp(1:250,:);
%[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_k_dec4)), 1:4096, 1024, hanning(127), 0);
%TFR1024_Bulle_4_han127_k = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_k = rtfr_tmp(1:250,:);
%[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_l_dec4)), 1:4096, 1024, hanning(127), 0);
%TFR1024_Bulle_4_han127_l = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_l = rtfr_tmp(1:250,:);


%TFR_Time = (0:length(Sig_Bulle_4_i_dec4)-2)/F_Ech*4;
%TFR_Freq1024 = (0:249)/1024*F_Ech/4;
%figure
%subplot(221);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_i); axis('xy'); caxis([1e-6 1e-5])
%title('Voie i TFTR1024 hanning127 : 110 kHz');
%subplot(222);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_j); axis('xy'); caxis([2e-6 2e-5])
%title('Voie j TFTR1024 hanning127 : 122 kHz');
%subplot(223);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_k); axis('xy'); caxis([1e-7 1e-6])
%title('Voie k TFTR1024 hanning127 : 122 kHz');
%subplot(224);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_l); axis('xy'); caxis([2e-6 2e-5])
%title('Voie l TFTR1024 hanning127 : 110 kHz');
%figure
%subplot(221);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_i); axis('xy'); caxis([1e-5 1e-4])
%title('Voie i RTFTR1024 hanning127 : 110 kHz');
%subplot(222);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_j); axis('xy'); caxis([1e-5 1e-4])
%title('Voie j RTFTR1024 hanning127 : 122 kHz');
%subplot(223);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_k); axis('xy'); caxis([1e-6 1e-5])
%title('Voie k RTFTR1024 hanning127 : 122 kHz');
%subplot(224);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_l); axis('xy'); caxis([1e-5 1e-4])
%title('Voie l RTFTR1024 hanning127 : 110 kHz');

% Extraction des Doppler
% TFRSP1024 
%Dop_TFR1024_Bulle_4_i=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_i(20:250,:))./sum(TFR1024_Bulle_4_han127_i(20:250,:));
%Dop_TFR1024_Bulle_4_j=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_j(20:250,:))./sum(TFR1024_Bulle_4_han127_j(20:250,:));
%Dop_TFR1024_Bulle_4_k=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_k(20:250,:))./sum(TFR1024_Bulle_4_han127_k(20:250,:));
%Dop_TFR1024_Bulle_4_l=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_l(20:250,:))./sum(TFR1024_Bulle_4_han127_l(20:250,:));
%figure
%subplot(221);plot(TFR_Time, Dop_TFR1024_Bulle_4_i); axis([0 0.25 0 max(TFR_Freq1024)]); grid
%title('Voie i TFTR1024 hanning127 : 110 kHz');
%subplot(222);plot(TFR_Time, Dop_TFR1024_Bulle_4_j); axis([0 0.25 0 max(TFR_Freq1024)]); grid
%title('Voie j TFTR1024 hanning127 : 122 kHz');
%subplot(223);plot(TFR_Time, Dop_TFR1024_Bulle_4_k); axis([0 0.25 0 max(TFR_Freq1024)]); grid
%title('Voie k TFTR1024 hanning127 : 122 kHz');
%subplot(224);plot(TFR_Time, Dop_TFR1024_Bulle_4_l); axis([0 0.25 0 max(TFR_Freq1024)]); grid
%title('Voie l TFTR1024 hanning127 : 110 kHz');
% RTFRSP1024
%Dop_RTFR1024_Bulle_4_i=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_i(20:250,:))./sum(RTFR1024_Bulle_4_han127_i(20:250,:));
%Dop_RTFR1024_Bulle_4_j=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_j(20:250,:))./sum(RTFR1024_Bulle_4_han127_j(20:250,:));
%Dop_RTFR1024_Bulle_4_k=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_k(20:250,:))./sum(RTFR1024_Bulle_4_han127_k(20:250,:));
%Dop_RTFR1024_Bulle_4_l=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_l(20:250,:))./sum(RTFR1024_Bulle_4_han127_l(20:250,:));

% MVA Mordant Pinton Michel
sig = conj(diff(detrend(Sig_Bulle_4_i_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.1;
wMVA;
Dop_MVA_Bulle_4_i = freq(1,:)*F_Ech/4;
InvHess_MVA_Bulle_4_i = iHes;


sig = conj(diff(detrend(Sig_Bulle_4_j_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.05;
wMVA;
Dop_MVA_Bulle_4_j = freq(1,:)*F_Ech/4;
InvHess_MVA_Bulle_4_j = iHes;

sig = conj(diff(detrend(Sig_Bulle_4_k_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.05;
wMVA;
Dop_MVA_Bulle_4_k = freq(1,:)*F_Ech/4;
InvHess_MVA_Bulle_4_k = iHes;

sig = conj(diff(detrend(Sig_Bulle_4_l_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
MVAinit=0.1;
wMVA;
Dop_MVA_Bulle_4_l = freq(1,:)*F_Ech/4;
InvHess_MVA_Bulle_4_l = iHes;
clear sig

REV_MVA_FLG = 0;
if (REV_MVA_FLG)
	sig = flipud(diff(detrend(Sig_Bulle_4_i_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_i = fliplr(freq(1,:))*F_Ech/4;
	InvHess_REV_MVA_Bulle_4_i = fliplr(iHes);
	
	sig = flipud(diff(detrend(Sig_Bulle_4_j_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_j = fliplr(freq(1,:))*F_Ech/4;
	InvHess_REV_MVA_Bulle_4_j = fliplr(iHes);
	
	sig = flipud(diff(detrend(Sig_Bulle_4_k_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_k = fliplr(freq(1,:))*F_Ech/4;
	InvHess_REV_MVA_Bulle_4_k = fliplr(iHes);
	
	sig = flipud(diff(detrend(Sig_Bulle_4_l_dec4((end+1)/2-1024-10:(end+1)/2+1024+11))));
	MVAinit=0.1;
	wMVA;
	Dop_REV_MVA_Bulle_4_l = fliplr(freq(1,:))*F_Ech/4;
	InvHess_REV_MVA_Bulle_4_l = fliplr(iHes);
end
figure
subplot(221); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_MVA_Bulle_4_i,'b'); axis([0 0.125 0 max(TFR_Freq1024)]); grid
subplot(222); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_MVA_Bulle_4_j,'g'); axis([0 0.125 0 max(TFR_Freq1024)]); grid
subplot(223); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_MVA_Bulle_4_k,'r'); axis([0 0.125 0 max(TFR_Freq1024)]); grid
subplot(224); plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_MVA_Bulle_4_l,'k'); axis([0 0.125 0 max(TFR_Freq1024)]); grid
% Overlay REV_MVA
subplot(221); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_REV_MVA_Bulle_4_i,'b-.'); hold off
subplot(222); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_REV_MVA_Bulle_4_j,'g-.'); hold off
subplot(223); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_REV_MVA_Bulle_4_k,'r-.'); hold off
subplot(224); hold on; plot((0:length(Dop_MVA_Bulle_4_i)-1)/F_Ech*4,Dop_REV_MVA_Bulle_4_l,'k-.'); hold off
% Overlay DOP_RTFR
I_Time_RTFR = 4096-1024:4096+1023;
subplot(221); hold on; plot((0:length(I_Time_RTFR)-1)/F_Ech*4,Dop_RTFR1024_Bulle_4_i(I_Time_RTFR),'b--'); hold off
subplot(222); hold on; plot((0:length(I_Time_RTFR)-1)/F_Ech*4,Dop_RTFR1024_Bulle_4_j(I_Time_RTFR),'g--'); hold off
subplot(223); hold on; plot((0:length(I_Time_RTFR)-1)/F_Ech*4,Dop_RTFR1024_Bulle_4_k(I_Time_RTFR),'r--'); hold off
subplot(224); hold on; plot((0:length(I_Time_RTFR)-1)/F_Ech*4,Dop_RTFR1024_Bulle_4_l(I_Time_RTFR),'k--'); hold off
subplot(221); title('Doppler (MVA REV-MVA RTFRSP Voie i Bulle #4')
subplot(222); title('Doppler (MVA REV-MVA RTFRSP Voie j Bulle #4')
subplot(223); title('Doppler (MVA REV-MVA RTFRSP Voie k Bulle #4')
subplot(224); title('Doppler (MVA REV-MVA RTFRSP Voie l Bulle #4')
% TFRCW sur 1024 points

tfr_tmp = tfrcw(detrend(conj(Sig_Bulle_4_i_dec4)), 1:4096, 1024, hanning(63), hanning(127), 1, 1);
TFRCW1024_Bulle_4_han127_i = tfr_tmp(1:500,:);
tfr_tmp = tfrcw(detrend(conj(Sig_Bulle_4_j_dec4)), 1:4096, 1024, hanning(63), hanning(127), 1, 1);
TFRCW1024_Bulle_4_han127_j = tfr_tmp(1:500,:);
tfr_tmp = tfrcw(detrend(conj(Sig_Bulle_4_k_dec4)), 1:4096, 1024, hanning(63), hanning(127), 1, 1);
TFRCW1024_Bulle_4_han127_k = tfr_tmp(1:500,:);
tfr_tmp = tfrcw(detrend(conj(Sig_Bulle_4_l_dec4)), 1:4096, 1024, hanning(63), hanning(127), 1, 1);
TFRCW1024_Bulle_4_han127_l = tfr_tmp(1:500,:);
