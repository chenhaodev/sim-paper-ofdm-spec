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
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
fs = 4;  % normlized sampling rate
fc = 2; %normlized RF-carrier rate
M = 20; % window length
iter = 100;
load signal_ofdm.mat
[s1, s2, s3] = size(signal_ofdm);
L = s1;

for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise1
		n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex, real + image
		% noise2
		x = n;
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
		%energy = abs(y).^2;  % feature-power
		%energy_fin(kk) =(1/L).*sum(energy); % test statistic (normlized)
        energy_fin(kk) = (1/length(L)).*norm(y);% feature-power, test statistic (normlized)
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
csd_thresh_est = thresh;
save ./Data/csd_thres.mat csd_thresh_est

