% calculate the threhsold in compressed sensing based cyclic feature detection by simulations. 
% 1. assume only noise is received, i.e., primary user is absent.
% false alarm, then
% 2. probability of false alarm = energy above threshold/No. of iteration.
% by chenhaomails@gmail.com
% ver 0.1

clc
close all
clear all
addpath('../Util/')
addpath('../Data/')
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
fs = 4;  % normlized sampling rate
fc = 2; %normlized RF-carrier rate
M = 10; % window length
iter = 100;
load signal_ofdm.mat
[s1, s2, s3] = size(signal_ofdm);
cs.ratio = 4;  % cs sampling ratio
L_origin = s1;
L = s1/cs.ratio;
cs.N = L_origin;
cs.M = round(cs.N/cs.ratio);
% ideal pass filter
lpf_temp = zeros(1,(cs.N-1)); lpf_temp(1) = 1;
lpf_z = interpft(lpf_temp,cs.N)/cs.N;
f_lpf = fft(lpf_z);

%% loop %%

for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise
		n=(randn(1,cs.N)+sqrt(-1)*randn(1,cs.N))./(sqrt(2)); % complex, real + image
		% input
		x_in = n;
        % cs sampling
		pn = pn_gen(cs.N); % generate pn sequence 
		mixed_sig = x_in.*pn;  % signal mixing
		f_sig = fft(mixed_sig);
		filter_out = ifft(f_lpf.*f_sig); % filtering via ideal pass filter
		cs.y = downsample(filter_out, cs.ratio); % low rate sampling
		x = cs.y; 
		% csd function (area power detection)
        [Spec, f, alpha] = cyclic_spectrum_new(x, cs.M, fs/cs.ratio, M, 'no-display');
		%feature range: f ~[-fs fs];
        f_mid = round(length(f)/2);
        f_index_range = 1:1:length(f);
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
cs_csd_thresh_est = thresh;
save ./cs_csd_thres.mat cs_csd_thresh_est

