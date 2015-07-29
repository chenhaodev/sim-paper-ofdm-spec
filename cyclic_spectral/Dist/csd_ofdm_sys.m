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
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
fs = 4;  % normlized sampling rate
fc = 2; %normlized RF-carrier rate
M = 20; % window length
iter = 100;
load ./Data/signal_ofdm.mat
L = length(signal_ofdm(:,1,1));
load ./Data/csd_thres.mat

for m = 1:length(Pf)
    m
    i = 0;
	for kk=1:iter % Number of Monte Carlo Simulations
		% noise
		n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex, real + image
		% pu signal
        s = signal_ofdm(:,1,kk)';
        x = sqrt(snr).*s + n; 
        % csd function (area power detection)
        [Spec, f, alpha] = cyclic_spectrum_new(x, L, fs, M, 'silent');
		%feature range: f ~[-fs/4 fs/4];
        f_mid = round(length(f)/2);
        f_index_range = (f_mid - floor(length(f)/(fs))+1) : (f_mid + floor(length(f)/(fs)));
        f_area = f(f_index_range);
        a_index_range = 1 : (round(length(alpha) * 0.5));
        a_area = alpha(a_index_range);
        % feature extraction
		y = (Spec(a_index_range,f_index_range));
        energy_fin = (1/length(L)).*norm(y);% feature-power, test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = csd_thresh_est(m); %use estimated threshold to detect ofdm
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
    end
    Pd_csd_ofdm_th_est(m) = i/kk; 
end
plot(Pf, Pd_csd_ofdm_th_est, '*r'); hold on; 
save ./Data/csd_ofdm_sys.mat Pd_csd_ofdm_th_est 
load ./Data/ed_ofdm_sys.mat
plot(Pf, Pd_rnd, 'ob'); hold on;
plot(Pf, Pd_ofdm, '.y'); hold on;
plot(Pf, Pd_ofdm_th_est, '.g'); hold on;
legend('ROC of FD for OFDM using estimated threshold','ROC of ED for random signal', 'ROC of ED for OFDM', 'ROC of ED for OFDM using estimated threshold');
xlabel('Pf');
ylabel('Pd');
