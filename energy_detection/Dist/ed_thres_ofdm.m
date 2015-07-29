% calculate the threhsold in energy detection by simulations. 
% 1. assume only noise is received, i.e., primary user is absent.
% false alarm, then
% 2. probability of false alarm = energy above threshold/No. of iteration.
% by chenhaomails@gmail.com
% ver 0.1

clc
close all
clear all
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
load ./Data/sig_ofdm/signal_ofdm.mat
L = size(signal_ofdm,1);
iter = 100;

for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex, real + image
		y = n; 
		% energy
		energy = abs(y).^2; 
		energy_fin(kk) =(1/L).*sum(energy); % test statistic (normlized)
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
thresh_est = thresh;
save ./Data/thres/threshold.mat thresh_est

