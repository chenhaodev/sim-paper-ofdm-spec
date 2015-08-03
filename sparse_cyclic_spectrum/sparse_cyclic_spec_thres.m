% calculate the threhsold in cyclic feature detection by simulations. 
% 1. assume only noise is received, i.e., primary user is absent.
% false alarm, then
% 2. probability of false alarm = energy above threshold/No. of iteration.
% by chenhaomails@gmail.com
% ver 0.1

% header

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

sig.type = 'fsk'; % 'fsk'
sig.fs = 1;
sig.M = 1;

if strcmpi(sig.type,'fsk') % default signal
	load fsk.mat             
else
	error('signal type not exist!!');
end

sig.x = fsk_real;
sig.N = length(sig.x);

snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR

Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
iter = 25; % Monte Carlo simulation

%% Loop %%
for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise1
		n=(randn(1,sig.N)+sqrt(-1)*randn(1,sig.N))./(sqrt(2)); % complex noise
		sig.x = n;
		% test-core
		[Spec_cs, f, alpha] = cyclic_xcorr(sig.x, sig.N, sig.fs, 'non-show'); 
		% feature-extract
		[feature_o, ~] = feature_test(Spec_cs, f, [-0.25 +0.25], alpha, [-0.15 +0.15]) 
		% energy 
        energy_fin(kk) = (1/length(sig.N)).*norm(feature_o);% test statistic (normlized)
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
cs_csd_thresh_est = thresh;
save ./Data/cs_csd_thres.mat cs_csd_thresh_est

