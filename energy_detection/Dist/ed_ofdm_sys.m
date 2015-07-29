% This code is to plot receiver operating characteristic curve for simple energy
% detection, when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is available
% analytically.
% Code written by: Sanket Kalamkar, Indian Institute of Technology Kanpur,
% India.


clc
close all
clear all
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
load ./Data/sig_ofdm/signal_ofdm.mat
L = size(signal_ofdm,1);
load ./Data/thres/threshold.mat 

%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf)
    m
    i = 0;
	ii = 0;
    iii = 0;
	for kk=1:100 % Number of Monte Carlo Simulations
		% signal generation
        %mask  = randi([0 1],1,ofdm.B);
		%if mask == 0
		%	dataIFFTGI = zeros(length(dataIFFTGI),1);
		%end
		n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex, real + image
        s = sqrt(snr).*(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2));
		%n = randn(1,L);
		%s = sqrt(snr).*randn(1,L); % pu signal
		y = s + n; 
		% energy collection
		energy = abs(y).^2; 
		energy_fin =(1/L).*sum(energy); % test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = (qfuncinv(Pf(m))./sqrt(L))+ 1; % threshold. It's determined by Pf, however, is not practical sometimes
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
		% ofdm signal energy detection
		s = sqrt(snr)./sqrt(2).*signal_ofdm(:,1,kk)'; % pu signal
		y = s + n; 
		% energy collection
		energy = abs(y).^2; 
		energy_fin =(1/L).*sum(energy); % test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = (qfuncinv(Pf(m))./sqrt(L))+ 1; % threshold. It's determined by Pf, however, is not practical sometimes
		if(energy_fin >= thresh(m))  
			ii = ii+1;
		end
		thresh(m) = thresh_est(m); %use estimated threshold to detect ofdm
		if(energy_fin >= thresh(m))  
			iii = iii+1;
		end
	end
Pd_rnd(m) = i/kk; 
Pd_ofdm(m) = ii/kk; 
Pd_ofdm_th_est(m) = iii/kk; 
end
plot(Pf, Pd_rnd, 'ob');
hold on
plot(Pf, Pd_ofdm, '*r');
hold on
plot(Pf, Pd_ofdm_th_est, '.g');
legend('ROC of random signal', 'ROC of ofdm signal', 'ROC of ofdm using estimated threshold');
save ed_ofdm_sys.mat Pd_rnd Pd_ofdm Pd_ofdm_th_est