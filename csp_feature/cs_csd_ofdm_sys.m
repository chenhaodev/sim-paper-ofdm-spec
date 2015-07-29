% This code is to plot receiver operating characteristic curve for compressed sensing based
% cyclic feature detection, when the primary signal is ofdm modulated signal noise is
% addive white real Gaussian. Here, the threshold is available analytically.
% author: chenhaomails@gmail.com
% 2015.7

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
% load trained threshold
load cs_csd_thres.mat

for m = 1:length(Pf)
    m
    i = 0;
	for kk=1:iter % Number of Monte Carlo Simulations
		% noise
		n=(randn(1,cs.N)+sqrt(-1)*randn(1,cs.N))./(sqrt(2)); % complex, real + image
		% pu signal
        s = signal_ofdm(:,1,kk)';
        x_in = sqrt(snr).*s + n; 
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
        energy_fin = (1/length(L)).*norm(y);% feature-power, test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = cs_csd_thresh_est(m); %use estimated threshold to detect ofdm
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
    end
    Pd_cs_csd_ofdm_th_est(m) = i/kk; 
end
plot(Pf, Pd_cs_csd_ofdm_th_est, '*r'); hold on; 
save ./csd_ofdm_sys.mat Pd_cs_csd_ofdm_th_est 
load ed_ofdm_sys.mat
plot(Pf, Pd_rnd, 'ob'); hold on;
plot(Pf, Pd_ofdm, '.y'); hold on;
plot(Pf, Pd_ofdm_th_est, '.g'); hold on;
legend('ROC of FD for OFDM using estimated threshold','ROC of ED for random signal', 'ROC of ED for OFDM', 'ROC of ED for OFDM using estimated threshold');
xlabel('Pf');
ylabel('Pd');
