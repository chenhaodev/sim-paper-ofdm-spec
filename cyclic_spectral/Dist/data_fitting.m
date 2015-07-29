
% cubic fitting

clc; clear; close all;

load ./Data/csd_ofdm_sys.mat
load ./Data/ed_ofdm_sys.mat

Pf = 0.01:0.01:1;

  p11 = 1.9842;
  p12 = -4.0262;
  p13 = 2.8791;
  p14 = 0.20618;

x1 = Pd_csd_ofdm_th_est;
y1 = p11*(x1.^3) + p12*(x1.^2) + p13*x1 + p14 ;
plot(x1, y1, '.b');
hold on;


  p21 = 1.486;
  p22 = -3.255;
  p23 = 2.5915;
  p24 = 0.20247;

x2 = Pd_ofdm_th_est;
y2 = p21*(x2.^3) + p22*(x2.^2) + p23*x2 + p24 ;
plot(x2, y2, '.r');
hold on;

legend('ROC of FD for OFDM using estimated threshold','ROC of ED for OFDM using estimated threshold');
xlabel('Pf');
ylabel('Pd');