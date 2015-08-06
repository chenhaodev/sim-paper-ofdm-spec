% This code is to plot receiver operating characteristic curve for
% cs based detection for cyclostationary signal,
% , when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is pre-estimated.
% author: chenhaomails@gmail.com
% update: 15/08/05

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

% Header 

sig.type = 'fsk'; % 'fsk'
sig.fs = 1;
sig.M = 1;

if strcmpi(sig.type,'fsk') % default signal
	load fsk.mat             
else
	error('signal type not exist!!');
end

sig.x = fsk_real(1:64);
sig.x = sig.x ./ norm(sig.x);
%sig.x = fsk_real(65:128);
sig.N = length(sig.x);

Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
iter = 20; % Monte Carlo simulation

Pf = 0.01:0.01:1; % Pf = Probability of False Alarm

snr_dB = 5; % SNR in decibels, NOTE !!!!
snr = 10.^(snr_dB./10); % Linear Value of SNR
load ./Data/thres_sparse_cyclic_spec.mat %thres_sparse_cyclic_spec_est

for m = 1:length(Pf)
    m
    i = 0;
	for kk=1:iter % Number of Monte Carlo Simulations
		% noise
		n= randn(1,sig.N); 
		% pu signal
        s = sig.x;
        x = sqrt(snr).*s + n; 
		% test, cs_cyc_spec + feature extract
		[hat_spec] = sparse_cyclic_spec(sig.x, sig.N, sig.fs, 'non-show');
		[out] = feature_extract(abs(hat_spec), 1:sig.N, 0.2, 1:sig.N, 0.2);
        energy_fin = (1/length(L)).*norm(out);% feature-power, test statistic (normlized)
		thresh(m) = thres_sparse_cyclic_spec_est(m); %use estimated threshold 
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
    end
    Pd_sparse_cyclic_spec(m) = i/kk; 
end
plot(Pf, Pd_sparse_cyclic_spec, '*r'); hold on; 
save ./Data/Pd_sparse_cyclic_spec.mat Pd_sparse_cyclic_spec 
% energy(object) is different, make them all detecting fsk signal, then compare them @ the same SNR
%load ./Data/ed_ofdm_sys.mat
%plot(Pf, Pd_rnd, 'ob'); hold on;
%plot(Pf, Pd_ofdm, '.y'); hold on;
%plot(Pf, Pd_ofdm_th_est, '.g'); hold on;
%legend('ROC of FD for OFDM using estimated threshold','ROC of ED for random signal', 'ROC of ED for OFDM', 'ROC of ED for OFDM using estimated threshold');
xlabel('Pf');
ylabel('Pd');
