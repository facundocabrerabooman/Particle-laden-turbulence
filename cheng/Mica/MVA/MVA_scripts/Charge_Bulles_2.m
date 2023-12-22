% Charge_Bulles.m

FILE_NAME = 'bulle_001.mac';
[devs, data]=read_mac(FILE_NAME);
F_Ech = devs{1}.sampling;
F_Cent = devs{1}.modulation;
L_Sig = length(data{1});
Cson = 340;
Psi = 26; Phi = atan(sin(pi/180*Psi))*180/pi;
Theta_Scatt_div2 = Phi;		% Theta = Pi - 2*Phi
% Parametres Acquisition
% Correction CB (8 juin 2005)
%data{1} : k autour de 110kHz, j autour de 122kHz
%data{2} : i autour de 110kHz, l autour de 122kHz

F_inc1 = 110000; F_inc2 = 122000;
F_inc_k = 110000;
F_inc_j = 122000;
F_inc_l = 122000;
F_inc_i = 110000;

Dop_i = Cson/(2*sin(pi/180*Phi))/F_inc_i;
Dop_j = Cson/(2*sin(pi/180*Phi))/F_inc_j;
Dop_k = Cson/(2*sin(pi/180*Phi))/F_inc_k;
Dop_l = Cson/(2*sin(pi/180*Phi))/F_inc_l;

% Demodulation + Spectres
NFFT = 2048; NOVERLAP = NFFT/2; REM_AVG = 1; REM_ALIAS = 1; DIFF_FLG = 0;
Voie_l = data{1}.*exp(i*2*pi*(0:L_Sig-1)'*(F_Cent-F_inc1)/F_Ech);
Voie_j = data{1}.*exp(i*2*pi*(0:L_Sig-1)'*(F_Cent-F_inc2)/F_Ech);
Voie_i = data{2}.*exp(i*2*pi*(0:L_Sig-1)'*(F_Cent-F_inc1)/F_Ech);
Voie_k = data{2}.*exp(i*2*pi*(0:L_Sig-1)'*(F_Cent-F_inc2)/F_Ech);
[P_i, F, N_AVG] = mexSpectrumC(Voie_i, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech);
P_j = mexSpectrumC(Voie_j, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech);
P_k = mexSpectrumC(Voie_k, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech);
P_l = mexSpectrumC(Voie_l, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech);
figure
subplot(221);plot(F, 10*log10(P_i), 'b');axis([-20000 20000 -85 -15]);grid;title('Voie i : F_{inc} = 110 kHz')
subplot(222);plot(F, 10*log10(P_j), 'g');axis([-20000 20000 -85 -15]);grid;title('Voie j : F_{inc} = 122 kHz')
subplot(223);plot(F, 10*log10(P_k), 'r');axis([-20000 20000 -85 -15]);grid;title('Voie k : F_{inc} = 122 kHz')
subplot(224);plot(F, 10*log10(P_l), 'k');axis([-20000 20000 -85 -15]);grid;title('Voie l : F_{inc} = 110 kHz')

% Decimation x2
DECIMx2_flg =0;
if (DECIMx2_flg)
	Voie_i_dec2 = decimate(Voie_i, 2);
	Voie_j_dec2 = decimate(Voie_j, 2);
	Voie_k_dec2 = decimate(Voie_k, 2);
	Voie_l_dec2 = decimate(Voie_l, 2);
	[P_i_dec2, F_dec2, N_AVG] = mexSpectrumC(Voie_i_dec2, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/2);
	P_j_dec2 = mexSpectrumC(Voie_j_dec2, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/2);
	P_k_dec2 = mexSpectrumC(Voie_k_dec2, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/2);
	P_l_dec2 = mexSpectrumC(Voie_l_dec2, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/2);
	figure
	subplot(221);plot(F_dec2, 10*log10(P_i_dec2), 'b');axis([-5000 5000 -85 -15]);grid;title('Voie i : F_{inc} = 110 kHz')
	subplot(222);plot(F_dec2, 10*log10(P_j_dec2), 'g');axis([-5000 5000 -85 -15]);grid;title('Voie j : F_{inc} = 122 kHz')
	subplot(223);plot(F_dec2, 10*log10(P_k_dec2), 'r');axis([-5000 5000 -85 -15]);grid;title('Voie k : F_{inc} = 122 kHz')
	subplot(224);plot(F_dec2, 10*log10(P_l_dec2), 'k');axis([-5000 5000 -85 -15]);grid;title('Voie l : F_{inc} = 110 kHz')
end

% Decimation x4
Voie_i_dec4 = decimate(Voie_i, 4);
Voie_j_dec4 = decimate(Voie_j, 4);
Voie_k_dec4 = decimate(Voie_k, 4);
Voie_l_dec4 = decimate(Voie_l, 4);
[P_i_dec4, F_dec4, N_AVG] = mexSpectrumC(Voie_i_dec4, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/4);
P_j_dec4 = mexSpectrumC(Voie_j_dec4, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/4);
P_k_dec4 = mexSpectrumC(Voie_k_dec4, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/4);
P_l_dec4 = mexSpectrumC(Voie_l_dec4, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/4);
figure
subplot(221);plot(F_dec4, 10*log10(P_i_dec4), 'b');axis([-5000 5000 -90 -30]);grid;title('Voie i Decimx4 : F_{inc} = 110 kHz')
subplot(222);plot(F_dec4, 10*log10(P_j_dec4), 'g');axis([-5000 5000 -90 -30]);grid;title('Voie j Decimx4 : F_{inc} = 122 kHz')
subplot(223);plot(F_dec4, 10*log10(P_k_dec4), 'r');axis([-5000 5000 -90 -30]);grid;title('Voie k Decimx4 : F_{inc} = 122 kHz')
subplot(224);plot(F_dec4, 10*log10(P_l_dec4), 'k');axis([-5000 5000 -90 -30]);grid;title('Voie l Decimx4 : F_{inc} = 110 kHz')

% Decimation x8
DECIMx8_flg =0;
if (DECIMx8_flg)
	Voie_i_dec8 = decimate(Voie_i, 8);
	Voie_j_dec8 = decimate(Voie_j, 8);
	Voie_k_dec8 = decimate(Voie_k, 8);
	Voie_l_dec8 = decimate(Voie_l, 8);
	[P_i_dec8, F_dec8, N_AVG] = mexSpectrumC(Voie_i_dec8, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/8);
	P_j_dec8 = mexSpectrumC(Voie_j_dec8, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/8);
	P_k_dec8 = mexSpectrumC(Voie_k_dec8, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/8);
	P_l_dec8 = mexSpectrumC(Voie_l_dec8, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech/8);
	figure
	subplot(221);plot(F_dec8, 10*log10(P_i_dec8), 'b');grid
	subplot(222);plot(F_dec8, 10*log10(P_j_dec8), 'g');grid
	subplot(223);plot(F_dec8, 10*log10(P_k_dec8), 'r');grid
	subplot(224);plot(F_dec8, 10*log10(P_l_dec8), 'k');grid
end

% Voir Bulles : decimatin x4

% All :
figure
subplot(411);plot((0:262000)/F_Ech*4,abs(diff(Voie_i_dec4(1:262000+2))).^2); grid
axis([0 12 0 2e-6]); title('Voie i Decimx4 : F_{inc} = 110 kHz')
hold on; 
plot([22750 22750]/F_Ech*4,[0 2e-6],'r--');
plot([35400 35400]/F_Ech*4,[0 2e-6],'r--');
plot([84000 84000]/F_Ech*4,[0 2e-6],'r--');
plot([103000 103000]/F_Ech*4,[0 2e-6],'r--');
plot([132600 132600]/F_Ech*4,[0 2e-6],'r--');
plot([186100 186100]/F_Ech*4,[0 2e-6],'r--');
hold off
subplot(412);plot((0:262000)/F_Ech*4,abs(diff(Voie_j_dec4(1:262000+2))).^2); grid
axis([0 12 0 2e-6]); title('Voie j Decimx4 : F_{inc} = 122 kHz')
hold on; 
plot([22750 22750]/F_Ech*4,[0 2e-6],'r--');
plot([35400 35400]/F_Ech*4,[0 2e-6],'r--');
plot([84000 84000]/F_Ech*4,[0 2e-6],'r--');
plot([103000 103000]/F_Ech*4,[0 2e-6],'r--');
plot([132600 132600]/F_Ech*4,[0 2e-6],'r--');
plot([186100 186100]/F_Ech*4,[0 2e-6],'r--');
hold off
subplot(413);plot((0:262000)/F_Ech*4,abs(diff(Voie_k_dec4(1:262000+2))).^2); grid
axis([0 12 0 2e-6]); title('Voie k Decimx4 : F_{inc} = 122 kHz')
hold on; 
plot([22750 22750]/F_Ech*4,[0 2e-6],'r--');
plot([35400 35400]/F_Ech*4,[0 2e-6],'r--');
plot([84000 84000]/F_Ech*4,[0 2e-6],'r--');
plot([103000 103000]/F_Ech*4,[0 2e-6],'r--');
plot([132600 132600]/F_Ech*4,[0 2e-6],'r--');
plot([186100 186100]/F_Ech*4,[0 2e-6],'r--');
hold off
subplot(414);plot((0:262000)/F_Ech*4,abs(diff(Voie_l_dec4(1:262000+2))).^2); grid
axis([0 12 0 2e-6]); title('Voie l Decimx4 : F_{inc} = 110 kHz')
hold on; 
plot([22750 22750]/F_Ech*4,[0 2e-6],'r--');
plot([35400 35400]/F_Ech*4,[0 2e-6],'r--');
plot([84000 84000]/F_Ech*4,[0 2e-6],'r--');
plot([103000 103000]/F_Ech*4,[0 2e-6],'r--');
plot([132600 132600]/F_Ech*4,[0 2e-6],'r--');
plot([186100 186100]/F_Ech*4,[0 2e-6],'r--');
hold off

% Bulle_1 @ 22750
figure
Bulle_1 = 22750; Pos = Bulle_1;
subplot(411);plot((0:4095)/F_Ech*4,abs(diff(Voie_i_dec4(Pos-2048:Pos+2048))).^2); grid
title(sprintf('Bulle #1 Voie i Decimx4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(412);plot((0:4095)/F_Ech*4,abs(diff(Voie_j_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(413);plot((0:4095)/F_Ech*4,abs(diff(Voie_k_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(414);plot((0:4095)/F_Ech*4,abs(diff(Voie_l_dec4(Pos-2048:Pos+2048))).^2); grid

% Bulle_2 @ 35400
figure
Bulle_2 = 35400; Pos = Bulle_2;
subplot(411);plot((0:4095)/F_Ech*4,abs(diff(Voie_i_dec4(Pos-2048:Pos+2048))).^2); grid
title(sprintf('Bulle #2 Voie i Decimx4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(412);plot((0:4095)/F_Ech*4,abs(diff(Voie_j_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(413);plot((0:4095)/F_Ech*4,abs(diff(Voie_k_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(414);plot((0:4095)/F_Ech*4,abs(diff(Voie_l_dec4(Pos-2048:Pos+2048))).^2); grid

% Bulle_3 @ 84000
figure
Bulle_3 = 84000; Pos = Bulle_3;
subplot(411);plot((0:4095)/F_Ech*4,abs(diff(Voie_i_dec4(Pos-2048:Pos+2048))).^2); grid
title(sprintf('Bulle #3 Voie i Decimx4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(412);plot((0:4095)/F_Ech*4,abs(diff(Voie_j_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(413);plot((0:4095)/F_Ech*4,abs(diff(Voie_k_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(414);plot((0:4095)/F_Ech*4,abs(diff(Voie_l_dec4(Pos-2048:Pos+2048))).^2); grid

% Bulle_4 @ 103000
figure
Bulle_4 = 103000; Pos = Bulle_4;
subplot(411);plot((0:4095)/F_Ech*4,abs(diff(Voie_i_dec4(Pos-2048:Pos+2048))).^2); grid
title(sprintf('Bulle #4 Voie i Decimx4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(412);plot((0:4095)/F_Ech*4,abs(diff(Voie_j_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(413);plot((0:4095)/F_Ech*4,abs(diff(Voie_k_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(414);plot((0:4095)/F_Ech*4,abs(diff(Voie_l_dec4(Pos-2048:Pos+2048))).^2); grid

% Bulle_5 @ 132600
figure
Bulle_5 = 132600; Pos = Bulle_5;
subplot(411);plot((0:4095)/F_Ech*4,abs(diff(Voie_i_dec4(Pos-2048:Pos+2048))).^2); grid
title(sprintf('Bulle #5 Voie i Decimx4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(412);plot((0:4095)/F_Ech*4,abs(diff(Voie_j_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(413);plot((0:4095)/F_Ech*4,abs(diff(Voie_k_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(414);plot((0:4095)/F_Ech*4,abs(diff(Voie_l_dec4(Pos-2048:Pos+2048))).^2); grid

% Bulle_6 @ 186100
figure
Bulle_6 = 186100; Pos = Bulle_6;
subplot(411);plot((0:4095)/F_Ech*4,abs(diff(Voie_i_dec4(Pos-2048:Pos+2048))).^2); grid
title(sprintf('Bulle #6 Voie i Decimx4 : Indice #%d, Date %2.3f s', Pos, (Pos-1)/F_Ech*4));
subplot(412);plot((0:4095)/F_Ech*4,abs(diff(Voie_j_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(413);plot((0:4095)/F_Ech*4,abs(diff(Voie_k_dec4(Pos-2048:Pos+2048))).^2); grid
subplot(414);plot((0:4095)/F_Ech*4,abs(diff(Voie_l_dec4(Pos-2048:Pos+2048))).^2); grid

% Spectres Bulles All
window = hanning(4096); KMU = 6*norm(window)^2; F_Bulles = (-2048:2047)'/4096*F_Ech/4;
% Voie i
Pos = Bulle_1;
P_bulles = abs(fftshift(fft(detrend(Voie_i_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_2;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_i_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_3;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_i_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_4;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_i_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_5;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_i_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_6;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_i_dec4(Pos-2048:Pos+2047)).*window))).^2;
P_bulles_i = P_bulles/KMU;
% voie j
Pos = Bulle_1;
P_bulles = abs(fftshift(fft(detrend(Voie_j_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_2;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_j_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_3;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_j_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_4;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_j_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_5;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_j_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_6;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_j_dec4(Pos-2048:Pos+2047)).*window))).^2;
P_bulles_j = P_bulles/KMU;
% voie k
Pos = Bulle_1;
P_bulles = abs(fftshift(fft(detrend(Voie_k_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_2;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_k_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_3;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_k_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_4;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_k_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_5;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_k_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_6;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_k_dec4(Pos-2048:Pos+2047)).*window))).^2;
P_bulles_k = P_bulles/KMU;
% voie l
Pos = Bulle_1;
P_bulles = abs(fftshift(fft(detrend(Voie_l_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_2;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_l_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_3;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_l_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_4;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_l_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_5;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_l_dec4(Pos-2048:Pos+2047)).*window))).^2;
Pos = Bulle_6;
P_bulles = P_bulles + abs(fftshift(fft(detrend(Voie_l_dec4(Pos-2048:Pos+2047)).*window))).^2;
P_bulles_l = P_bulles/KMU;

figure
subplot(221);plot(F_Bulles, 10*log10(P_bulles_i),'b');grid;title('FFT 4096 Bulles 1-6 Voie i Decimx4 : 110 kHz'); axis([-5000 5000 -90 -30]);
subplot(222);plot(F_Bulles, 10*log10(P_bulles_j),'g');grid;title('FFT 4096 Bulles 1-6 Voie j Decimx4 : 122 kHz'); axis([-5000 5000 -90 -30]);
subplot(223);plot(F_Bulles, 10*log10(P_bulles_k),'r');grid;title('FFT 4096 Bulles 1-6 Voie k Decimx4 : 122 kHz'); axis([-5000 5000 -90 -30]);
subplot(224);plot(F_Bulles, 10*log10(P_bulles_l),'k');grid;title('FFT 4096 Bulles 1-6 Voie l Decimx4 : 110 kHz'); axis([-5000 5000 -90 -30]);

% Spectres No Bulles
NFFT = 2048; NOVERLAP = NFFT/2; window = hanning(NFFT);
% Voie i
X = Voie_i_dec4(1:Bulle_1-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_AVG; P_no_bulles_i = P_tmp*N_AVG;
X = Voie_i_dec4(Bulle_1+2048:Bulle_2-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_i = P_no_bulles_i + P_tmp*N_AVG;
X = Voie_i_dec4(Bulle_2+2048:Bulle_3-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_i = P_no_bulles_i + P_tmp*N_AVG;
X = Voie_i_dec4(Bulle_3+2048:Bulle_4-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_i = P_no_bulles_i + P_tmp*N_AVG;
X = Voie_i_dec4(Bulle_4+2048:Bulle_5-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_i = P_no_bulles_i + P_tmp*N_AVG;
X = Voie_i_dec4(Bulle_5+2048:Bulle_6-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_i = P_no_bulles_i + P_tmp*N_AVG;
X = Voie_i_dec4(Bulle_6+2048:end);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_i = P_no_bulles_i + P_tmp*N_AVG;
P_no_bulles_i = P_no_bulles_i/N_TOT;

% voie j
X = Voie_j_dec4(1:Bulle_1-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_AVG; P_no_bulles_j = P_tmp*N_AVG;
X = Voie_j_dec4(Bulle_1+2048:Bulle_2-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_j = P_no_bulles_j + P_tmp*N_AVG;
X = Voie_j_dec4(Bulle_2+2048:Bulle_3-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_j = P_no_bulles_j + P_tmp*N_AVG;
X = Voie_j_dec4(Bulle_3+2048:Bulle_4-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_j = P_no_bulles_j + P_tmp*N_AVG;
X = Voie_j_dec4(Bulle_4+2048:Bulle_5-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_j = P_no_bulles_j + P_tmp*N_AVG;
X = Voie_j_dec4(Bulle_5+2048:Bulle_6-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_j = P_no_bulles_j + P_tmp*N_AVG;
X = Voie_j_dec4(Bulle_6+2048:end);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_j = P_no_bulles_j + P_tmp*N_AVG;
P_no_bulles_j = P_no_bulles_j/N_TOT;

% voie k
X = Voie_k_dec4(1:Bulle_1-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_AVG; P_no_bulles_k = P_tmp*N_AVG;
X = Voie_k_dec4(Bulle_1+2048:Bulle_2-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_k = P_no_bulles_k + P_tmp*N_AVG;
X = Voie_k_dec4(Bulle_2+2048:Bulle_3-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_k = P_no_bulles_k + P_tmp*N_AVG;
X = Voie_k_dec4(Bulle_3+2048:Bulle_4-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_k = P_no_bulles_k + P_tmp*N_AVG;
X = Voie_k_dec4(Bulle_4+2048:Bulle_5-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_k = P_no_bulles_k + P_tmp*N_AVG;
X = Voie_k_dec4(Bulle_5+2048:Bulle_6-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_k = P_no_bulles_k + P_tmp*N_AVG;
X = Voie_k_dec4(Bulle_6+2048:end);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_k = P_no_bulles_k + P_tmp/N_AVG;
P_no_bulles_k = P_no_bulles_k/N_TOT;
% voie l
X = Voie_l_dec4(1:Bulle_1-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_AVG; P_no_bulles_l = P_tmp*N_AVG;
X = Voie_l_dec4(Bulle_1+2048:Bulle_2-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_l = P_no_bulles_l + P_tmp*N_AVG;
X = Voie_l_dec4(Bulle_2+2048:Bulle_3-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_l = P_no_bulles_l + P_tmp*N_AVG;
X = Voie_l_dec4(Bulle_3+2048:Bulle_4-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_l = P_no_bulles_l + P_tmp*N_AVG;
X = Voie_l_dec4(Bulle_4+2048:Bulle_5-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_l = P_no_bulles_l + P_tmp*N_AVG;
X = Voie_l_dec4(Bulle_5+2048:Bulle_6-2048);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_l = P_no_bulles_l + P_tmp*N_AVG;
X = Voie_l_dec4(Bulle_6+2048:end);
[P_tmp, F_no_bulles, N_AVG] = mexSpectrumC(X, NFFT, NOVERLAP, 1, 1, 0, window, F_Ech/4);
N_TOT = N_TOT+N_AVG; P_no_bulles_l = P_no_bulles_l + P_tmp*N_AVG;
P_no_bulles_l = P_no_bulles_l/N_TOT;

figure
subplot(221);plot(F_no_bulles, 10*log10(P_no_bulles_i),'b');grid;title('FFT 4096 No Bulles Voie i Decimx4 : 110 kHz'); axis([-5000 5000 -90 -30]);
subplot(222);plot(F_no_bulles, 10*log10(P_no_bulles_j),'g');grid;title('FFT 4096 No Bulles Voie j Decimx4 : 122 kHz'); axis([-5000 5000 -90 -30]);
subplot(223);plot(F_no_bulles, 10*log10(P_no_bulles_k),'r');grid;title('FFT 4096 No Bulles Voie k Decimx4 : 122 kHz'); axis([-5000 5000 -90 -30]);
subplot(224);plot(F_no_bulles, 10*log10(P_no_bulles_l),'k');grid;title('FFT 4096 No Bulles Voie l Decimx4 : 110 kHz'); axis([-5000 5000 -90 -30]);



% Test Temps-Frequence
% Bulle 4, Decimate x4
Bulle_4 = 103000; Pos = Bulle_4; s_Bulle_4_dec4 = Pos-2048:Pos+2048;
Sig_Bulle_4_i_dec4 = Voie_i_dec4(s_Bulle_4_dec4);
Sig_Bulle_4_j_dec4 = Voie_j_dec4(s_Bulle_4_dec4);
Sig_Bulle_4_k_dec4 = Voie_k_dec4(s_Bulle_4_dec4);
Sig_Bulle_4_l_dec4 = Voie_l_dec4(s_Bulle_4_dec4);
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

[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_i_dec4)), 1:4096, 1024, hanning(127), 0);
TFR1024_Bulle_4_han127_i = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_i = rtfr_tmp(1:250,:); 
[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_j_dec4)), 1:4096, 1024, hanning(127), 0);
TFR1024_Bulle_4_han127_j = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_j = rtfr_tmp(1:250,:);
[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_k_dec4)), 1:4096, 1024, hanning(127), 0);
TFR1024_Bulle_4_han127_k = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_k = rtfr_tmp(1:250,:);
[tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_l_dec4)), 1:4096, 1024, hanning(127), 0);
TFR1024_Bulle_4_han127_l = tfr_tmp(1:250,:); RTFR1024_Bulle_4_han127_l = rtfr_tmp(1:250,:);


TFR_Time = (0:length(Sig_Bulle_4_i_dec4)-2)/F_Ech*4;
TFR_Freq1024 = (0:249)/1024*F_Ech/4;
figure
subplot(221);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_i); axis('xy'); caxis([1e-6 1e-5])
title('Voie i TFTR1024 hanning127 : 110 kHz');
subplot(222);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_j); axis('xy'); caxis([2e-6 2e-5])
title('Voie j TFTR1024 hanning127 : 122 kHz');
subplot(223);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_k); axis('xy'); caxis([1e-7 1e-6])
title('Voie k TFTR1024 hanning127 : 122 kHz');
subplot(224);imagesc(TFR_Time, TFR_Freq1024, TFR1024_Bulle_4_han127_l); axis('xy'); caxis([2e-6 2e-5])
title('Voie l TFTR1024 hanning127 : 110 kHz');
figure
subplot(221);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_i); axis('xy'); caxis([1e-5 1e-4])
title('Voie i RTFTR1024 hanning127 : 110 kHz');
subplot(222);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_j); axis('xy'); caxis([1e-5 1e-4])
title('Voie j RTFTR1024 hanning127 : 122 kHz');
subplot(223);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_k); axis('xy'); caxis([1e-6 1e-5])
title('Voie k RTFTR1024 hanning127 : 122 kHz');
subplot(224);imagesc(TFR_Time, TFR_Freq1024, RTFR1024_Bulle_4_han127_l); axis('xy'); caxis([1e-5 1e-4])
title('Voie l RTFTR1024 hanning127 : 110 kHz');

% Extraction des Doppler
% TFRSP1024 
Dop_TFR1024_Bulle_4_i=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_i(20:250,:))./sum(TFR1024_Bulle_4_han127_i(20:250,:));
Dop_TFR1024_Bulle_4_j=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_j(20:250,:))./sum(TFR1024_Bulle_4_han127_j(20:250,:));
Dop_TFR1024_Bulle_4_k=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_k(20:250,:))./sum(TFR1024_Bulle_4_han127_k(20:250,:));
Dop_TFR1024_Bulle_4_l=(TFR_Freq1024(20:end)*TFR1024_Bulle_4_han127_l(20:250,:))./sum(TFR1024_Bulle_4_han127_l(20:250,:));
figure
subplot(221);plot(TFR_Time, Dop_TFR1024_Bulle_4_i); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie i TFTR1024 hanning127 : 110 kHz');
subplot(222);plot(TFR_Time, Dop_TFR1024_Bulle_4_j); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie j TFTR1024 hanning127 : 122 kHz');
subplot(223);plot(TFR_Time, Dop_TFR1024_Bulle_4_k); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie k TFTR1024 hanning127 : 122 kHz');
subplot(224);plot(TFR_Time, Dop_TFR1024_Bulle_4_l); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie l TFTR1024 hanning127 : 110 kHz');
% RTFRSP1024
Dop_RTFR1024_Bulle_4_i=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_i(20:250,:))./sum(RTFR1024_Bulle_4_han127_i(20:250,:));
Dop_RTFR1024_Bulle_4_j=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_j(20:250,:))./sum(RTFR1024_Bulle_4_han127_j(20:250,:));
Dop_RTFR1024_Bulle_4_k=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_k(20:250,:))./sum(RTFR1024_Bulle_4_han127_k(20:250,:));
Dop_RTFR1024_Bulle_4_l=(TFR_Freq1024(20:end)*RTFR1024_Bulle_4_han127_l(20:250,:))./sum(RTFR1024_Bulle_4_han127_l(20:250,:));
figure
subplot(221);plot(TFR_Time, Dop_RTFR1024_Bulle_4_i); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie i TFTR1024 hanning127 : 110 kHz');
subplot(222);plot(TFR_Time, Dop_RTFR1024_Bulle_4_j); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie j TFTR1024 hanning127 : 122 kHz');
subplot(223);plot(TFR_Time, Dop_RTFR1024_Bulle_4_k); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie k TFTR1024 hanning127 : 122 kHz');
subplot(224);plot(TFR_Time, Dop_RTFR1024_Bulle_4_l); axis([0 0.25 0 max(TFR_Freq1024)]); grid
title('Voie l TFTR1024 hanning127 : 110 kHz');
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


TFR_Time = (0:length(Sig_Bulle_4_i_dec4)-2)/F_Ech*4;
TFR_Freq1024 = (0:499)/1024*F_Ech/4/2;
figure
subplot(221);imagesc(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_i); axis('xy'); caxis([1e-6 1e-5])
title('Voie i TFRCW1024 hanning 63-127 sigma 1 : 110 kHz');
subplot(222);imagesc(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_j); axis('xy'); caxis([2.5e-6 2.5e-5])
title('Voie j TFRCW1024 hanning 63-127 sigma 1 : 122 kHz');
subplot(223);imagesc(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_k); axis('xy'); caxis([1e-7 1e-6])
title('Voie k TFRCW1024 hanning 63-127 sigma 1 : 122 kHz');
subplot(224);imagesc(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_l); axis('xy'); caxis([2e-6 2e-5])
title('Voie l TFRCW1024 hanning 63-127 sigma 1 : 110 kHz');
% figure
% subplot(221);pcolor(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_i); shading('flat'); caxis([1e-6 1e-5])
% title('Voie i TFRCW1024 hanning 63-127 sigma 1 : 110 kHz');
% subplot(222);pcolor(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_j); shading('flat'); caxis([2.5e-6 2.5e-5])
% title('Voie j TFRCW1024 hanning 63-127 sigma 1 : 122 kHz');
% subplot(223);pcolor(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_k); shading('flat'); caxis([1e-7 1e-6])
% title('Voie k TFRCW1024 hanning 63-127 sigma 1 : 122 kHz');
% subplot(224);pcolor(TFR_Time, TFR_Freq1024, TFRCW1024_Bulle_4_han127_l); shading('flat'); caxis([2e-6 2e-5])
% title('Voie l TFRCW1024 hanning 63-127 sigma 1 : 110 kHz');

% TFRSP sur 2048 points : out_of_memory !
%
% [tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_i_dec4)), 1:4096, 2048, hanning(127), 0);
% RTFR2048_Bulle_4_han127_i = rtfr_tmp(1:500,:);
% TFR2048_Bulle_4_han127_i = tfr_tmp(1:500,:);
% [tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_j_dec4)), 1:4096, 2048, hanning(127), 0);
% TFR2048_Bulle_4_han127_j = tfr_tmp(1:500,:);
% RTFR2048_Bulle_4_han127_j = rtfr_tmp(1:500,:);
% [tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_k_dec4)), 1:4096, 2048, hanning(127), 0);
% TFR2048_Bulle_4_han127_k = tfr_tmp(1:500,:);
% RTFR2048_Bulle_4_han127_k = rtfr_tmp(1:500,:);
% [tfr_tmp, rtfr_tmp] = mex_tfrrsp(detrend(conj(Sig_Bulle_4_l_dec4)), 1:4096, 2048, hanning(127), 0);
% TFR2048_Bulle_4_han127_l = tfr_tmp(1:500,:);
% RTFR2048_Bulle_4_han127_l = rtfr_tmp(1:500,:);
% TFR_Time = (0:length(Sig_Bulle_4_i_dec4)-2)/F_Ech*4;
% TFR_Freq2048 = (0:499)/2048*F_Ech/4;
% figure
% subplot(221);imagesc(TFR_Time, TFR_Freq2048, TFR2048_Bulle_4_han127_i);axis('xy')
% title('Voie i TFTR2048 hanning127 : 110 kHz');
% subplot(222);imagesc(TFR_Time, TFR_Freq2048, TFR2048_Bulle_4_han127_j);;axis('xy')
% title('Voie j TFTR2048 hanning127 : 122 kHz');
% subplot(223);imagesc(TFR_Time, TFR_Freq2048, TFR2048_Bulle_4_han127_k);;axis('xy')
% title('Voie k TFTR2048 hanning127 : 122 kHz');
% subplot(224);imagesc(TFR_Time, TFR_Freq2048, TFR2048_Bulle_4_han127_l);;axis('xy')
% title('Voie l TFTR2048 hanning127 : 110 kHz');
% figure
% subplot(221);imagesc(TFR_Time, TFR_Freq2048, RTFR2048_Bulle_4_han127_i);caxis([1e-5 1e-3]); axis('xy')
% title('Voie i RTFTR2048 hanning127 : 110 kHz');
% subplot(222);imagesc(TFR_Time, TFR_Freq2048, RTFR2048_Bulle_4_han127_j);caxis([1e-5 1e-3]); axis('xy')
% title('Voie j RTFTR2048 hanning127 : 122 kHz');
% subplot(223);imagesc(TFR_Time, TFR_Freq2048, RTFR2048_Bulle_4_han127_k);caxis([1e-6 1e-4]); axis('xy')
% title('Voie k RTFTR2048 hanning127 : 122 kHz');
% subplot(224);imagesc(TFR_Time, TFR_Freq2048, RTFR2048_Bulle_4_han127_l);caxis([1e-5 1e-3]); axis('xy')
% title('Voie l RTFTR2048 hanning127 : 110 kHz');

% Extraction des 3 composantes xyz de vitesse
Psi = 26*pi/180; Phi = atan(sin(Psi); Cson = 340;
F_inc_i = 110000; F_inc_j = 122000; F_inc_k = 110000; F_inc_l = 122000;

Dop2Vit_i = Cson/(2*sin(Phi))/F_inc_i;
Dop2Vit_j = Cson/(2*sin(Phi))/F_inc_j;
Dop2Vit_k = Cson/(2*sin(Phi))/F_inc_k;
Dop2Vit_l = Cson/(2*sin(Phi))/F_inc_l;

Vx =   (Dop2Vit_l*Dop_MVA_Bulle_4_l - Dop2Vit_k*Dop_MVA_Bulle_4_k)/(2*sin(Psi));
Vy =   (Dop2Vit_i*Dop_MVA_Bulle_4_i - Dop2Vit_j*Dop_MVA_Bulle_4_j)/(2*sin(Psi));
Vz_1 = (Dop2Vit_l*Dop_MVA_Bulle_4_l + Dop2Vit_k*Dop_MVA_Bulle_4_k)/(2*cos(Psi));
Vz_2 = (Dop2Vit_i*Dop_MVA_Bulle_4_i + Dop2Vit_j*Dop_MVA_Bulle_4_j)/(2*cos(Psi));


% Filtrages
Fc = (F_inc2 - F_inc1)/F_Ech*2;
[B,A] = butter(4, Fc);
help decimate
Diff_i = filter(B, A, my_diff(Voie_i));
P_Diff_i = mexSpectrumC(Voie_i, NFFT, NOVERLAP, REM_AVG, REM_ALIAS, DIFF_FLG, hanning(NFFT), F_Ech);
Sig_filt_mex=mexFilterC(B, A, F_trans/F_ech, Sig);
