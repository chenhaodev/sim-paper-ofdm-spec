% calculate the threhsold in cyclic feature detection by simulations. 
% 1. assume only noise is received, i.e., primary user is absent.
% false alarm, then
% 2. probability of false alarm = energy above threshold/No. of iteration.
% by chenhaomails@gmail.com
% ver 0.1

clc
close all
clear all


addpath('./Util/')
addpath('./Data/')

sig.type = 'fsk'; % 'default', 'rand', 'ofdm'
sig.cyclic = 16; 
sig.fs = 4; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 5; %window length

snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR

Pf = 0.01:0.01:1; % Pf = Probability of False Alarm

iter = 50;
load FSK_real.mat
x = fsk_real(1:512);
L = length(x);
%% Loop %%
for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise1
		n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex, real + image
		% noise2
		x = n;
        % thres est
        energy_fin(kk) = (1/length(L)).*norm(x);% test statistic (normlized)
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
ed_thres_est = thresh;
save ./Data/ed_thres_est.mat ed_thres_est

