%%
%When the rig falls we have one transitory g to 0g (1). After, we have a
%transitory 0g to deceleration (2).

clear all;clc; close all

Fs = 2990;

folderout = '/Volumes/landau1/TrCer_1000/test_transitory/';
cd(folderout)

% DDT particles
load('/Volumes/landau1/TrCer_1000/ddt/tracklongP_conc.mat')
tracklong_ddt = tracklongP_conc; clear tracklongP_conc

% % Dec particles
% load('/Volumes/landau1/TrCer_1000/dec/tracklongP_conc.mat')
% tracklong_dec = tracklongP_conc; clear tracklongP_conc

%% Order data

t_v_conc(:,1) = vertcat(tracklong_ddt.Tf);
t_v_conc(:,2) = vertcat(tracklong_ddt.Vy);
t_v_conc_sorted_ddt = sortrows(t_v_conc, 1); 
t_v_conc_ddt.Vy = t_v_conc_sorted_ddt(:,2);
t_v_conc_ddt.Tf = t_v_conc_sorted_ddt(:,1); clear t_v_conc

%% Plot PDFS first transitory

t1 = 0.4;

figure(1);hold on; clf

I = find(t_v_conc_ddt.Tf<t1);
data_t1_ddt_timefiltered.Tf = t_v_conc_ddt.Tf(I);
data_t1_ddt_timefiltered.Vy = t_v_conc_ddt.Vy(I);

pdfV_t1 = mkpdf5(data_t1_ddt_timefiltered,'Vy',100,10);

%%% plot particle vel
h1=semilogy(pdfV_t1.xpdf,pdfV_t1.pdf,'r-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfV_t1.mean,'r',LineWidth=3)

%%% Plot PDFS last transition

t2 = 2.1-0.4;

figure(1);hold on

I = find(t_v_conc_ddt.Tf>t2);
data_ts_ddt_timefiltered.Tf = t_v_conc_ddt.Tf(I);
data_ts_ddt_timefiltered.Vy = t_v_conc_ddt.Vy(I);

pdfV_ts = mkpdf5(data_ts_ddt_timefiltered,'Vy',100,10);

%%% plot particle vel
h3=semilogy(pdfV_ts.xpdf,pdfV_ts.pdf,'b-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfV_ts.mean,'b',LineWidth=3)

%%% Plot PDFS steady state

ts = [t1 t2];

figure(1);hold on

I = find(t_v_conc_ddt.Tf>ts(1) & t_v_conc_ddt.Tf<ts(2));
data_ts_ddt_timefiltered.Tf = t_v_conc_ddt.Tf(I);
data_ts_ddt_timefiltered.Vy = t_v_conc_ddt.Vy(I);

pdfV_ts = mkpdf5(data_ts_ddt_timefiltered,'Vy',100,10);

%%% plot particle vel
h2=semilogy(pdfV_ts.xpdf,pdfV_ts.pdf,'g-o',MarkerSize=5,LineWidth=2);hold on
xline(pdfV_ts.mean,'g',LineWidth=3)

%% 

xlabel('Vertical Velocity (mm/s)')
ylabel('PDF')
legend([h1 h2 h3],{'Data from 0s to 0.4s','0.4s to 1.7s','1.7s to 2.1s'})
xticks([-800:100:800])

savefig_FC('hist_transitory',8,6,'pdf')
savefig_FC('hist_transitory',8,6,'fig')









