% This code is to plot receiver operating characteristic curve for cyclic feature
% detection, when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is available analytically.
% author: chenhaomails@gmail.com
% 2015.7

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

% signal
load FSK_real.mat
x = fsk_real(1:512);

% threshold
load cs_csd_thres.mat

% attr for cs sampling and sparse reconstruction 
cs.sparse = 8;
cs.ratio = 8;
cs.iter = 32;
cs.N = length(x);
cs.M = round(cs.N/cs.ratio);
Phi = randn(cs.M,cs.N);

%% Loop %%

for m = 1:length(Pf)
    m
    i = 0;
	for kk=1:iter % Number of Monte Carlo Simulations
		% noise
		n=(randn(1,cs.N)+sqrt(-1)*randn(1,cs.N))./(sqrt(2)); % complex, real + image
		% pu signal
        s = x;
        x = sqrt(snr).*s + n; 
        % cs sampling
        y = Phi*x';
        A = Phi*(dftmtx(cs.N))^-1;
        [recov,~] = cosamp(y, A, cs.sparse, cs.iter);  
        hatx = (dftmtx(cs.N))^-1 * recov;
		% cs based csd function
		[Spec, f, alpha] = my_cyclic_spectrum(hatx', cs.N, sig.fs, sig.M,'no-show', 'cs_based_cyclic_spectrum');
		% main component
		Spec_r = reshape(Spec, 1, length(f)* length(alpha));
		sorted_spec = sort(Spec_r, 'descend'); 
		main_spec = sorted_spec(1:10);
		rest_spec = sorted_spec(11:end);
		%test_stat = norm(main_spec) ./ norm(rest_spec);
        energy_fin = (1/length(cs.N)).*norm(main_spec);% test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = cs_csd_thresh_est(m); %use estimated threshold to detect ofdm
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
    end
    Pd_cs_fsk_th_est(m) = i/kk; 
end
plot(Pf, Pd_cs_fsk_th_est, '*r'); hold on; 
save ./Data/cs_fsk_detect_sys.mat Pd_cs_fsk_th_est 
load ./Data/ed_ofdm_sys.mat
plot(Pf, Pd_rnd, 'ob'); hold on;
plot(Pf, Pd_ofdm, '.y'); hold on;
plot(Pf, Pd_ofdm_th_est, '.g'); hold on;
legend('ROC of CS-FD for FSK using estimated threshold','ROC of ED for random signal', 'ROC of ED for OFDM', 'ROC of ED for OFDM using estimated threshold');
xlabel('Pf');
ylabel('Pd');
