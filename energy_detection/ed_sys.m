% This code is to plot receiver operating characteristic curve for simple energy
% detection, when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is available
% analytically.
% Code written by: Sanket Kalamkar, Indian Institute of Technology Kanpur,
% India.


clc
close all
clear all
L = 100;
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf)
    m
    i = 0;
	for kk=1:100 % Number of Monte Carlo Simulations
		% signal generation
		n = randn(1,L);
		s = sqrt(snr).*randn(1,L); % pu signal
		y = s + n; 
		% energy collection
		energy = abs(y).^2; 
		energy_fin =(1/L).*sum(energy); % test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = (qfuncinv(Pf(m))./sqrt(L))+ 1; % threshold. It's determined by Pf, however, is not practical sometimes
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
	end
Pd(m) = i/kk; 
end
plot(Pf, Pd)
hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.
thresh = (qfuncinv(Pf)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
plot(Pf, Pd_the, 'r')
hold on
