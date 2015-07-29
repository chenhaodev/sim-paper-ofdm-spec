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
cs.ratio = 10;
cs.N = L;
cs.M = round(cs.N/cs.ratio);

for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise1
		n= (randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex, real + image
		% noise2
		x = n;
		% cs sampling 
		pn = pn_gen(cs.N);% signal mixing
		mixed_sig = x.*pn;  
		lpf_temp = zeros(1,(cs.N-1)); % bandpass: side-band
		lpf_temp(1) = 1; % ideal pass filter
		lpf_z = interpft(lpf_temp,cs.N)/cs.N; % impulse response
		f_lpf = fft(lpf_z);
		f_sig = fft(mixed_sig);
		filter_out = ifft(f_lpf.*f_sig); %multiplication in freq, = convolution in time
		cs.y = downsample(filter_out, cs.ratio); %defactor = downsample ratio, low pass filtering and downsampling
		%feature analysis (acf_constell)
		% 1. extract/set apart the significant point from autocorrelated constellation
		% 2. calculate the mean location of the left points (non-significant point,mean_g)
		% 3. feature = distance(mean_g, significant point)
		o1 = mean(real(cs.y));
		o2 = mean(imag(cs.y));
		for mm = 1:1:length(cs.y)
		    temp = cs.y(mm); 
		    a1 = real(temp);
		    a2 = imag(temp);
		    dist_point(mm) = sqrt((a1-o1)^2 + (a2-o2)^2);
		end
		[dist_val,dis_inx] = sort(dist_point,'descend');
		key_inx = dis_inx(1);  
		group_inx = setxor(dis_inx, dis_inx(1));
		mean_g_r = mean(real(cs.y(group_inx))); % non-significant point
		mean_g_i = mean(imag(cs.y(group_inx)));
		key_g_r = real(cs.y(key_inx)); % significant point
		key_g_i = imag(cs.y(key_inx));
		acf_constell = sqrt((key_g_r-mean_g_r)^2 + (key_g_i-mean_g_i)^2);  %feature

		%energy = abs(y).^2;  % feature-power
		%energy_fin(kk) =(1/L).*sum(energy); % test statistic (normlized)
        energy_fin(kk) = (1/length(L)).*acf_constell;% feature-power, test statistic (normlized)
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
acf_constell_thresh_est = thresh;
save ./Data/acf_constell_thres.mat acf_constell_thresh_est

