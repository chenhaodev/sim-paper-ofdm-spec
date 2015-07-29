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
L = length(x);

% threshold
load ed_thres_est.mat

%% Loop %%

for m = 1:length(Pf)
    m
    i = 0;
	for kk=1:iter % Number of Monte Carlo Simulations
		% noise
		n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex, real + image
		% pu signal
        s = x;
        x = sqrt(snr).*s + n; 
        % energy detect
        energy_fin = (1/length(L)).*norm(x);% test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = ed_thres_est(m); %use estimated threshold to detect ofdm
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
    end
    Pd_ed_fsk_th_est(m) = i/kk; 
end
plot(Pf, Pd_ed_fsk_th_est, 'ob'); hold on; 
save ./Data/Pd_ed_fsk_th_est.mat
load ./Data/cs_fsk_detect_sys.mat  
plot(Pf, Pd_cs_fsk_th_est, '*r'); hold on;
legend('ROC of ED for FSK signal' ,'ROC of CS-FD for FSK using estimated threshold');
xlabel('Pf');
ylabel('Pd');
