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

sig.type = 'fsk'; % 'default', 'rand', 'ofdm'
sig.cyclic = 16; 
sig.fs = 4; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 5; %window length

snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR

Pf = 0.01:0.01:1; % Pf = Probability of False Alarm

iter = 100;
load FSK_real.mat
x = fsk_real(1:512);

% attr for cs sampling and sparse reconstruction 
cs.sparse = 16;
cs.ratio = 4;
cs.iter = 10;
cs.N = length(x);
cs.M = round(cs.N/cs.ratio);
Phi = randn(cs.M,cs.N);

%% Loop %%
for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise1
		n=(randn(1,cs.N)+sqrt(-1)*randn(1,cs.N))./(sqrt(2)); % complex, real + image
		% noise2
		x = n;
        % cs sampling
        y = Phi*x';
        A = Phi*(dctmtx(cs.N))^-1;
        [recov,~] = cosamp(y, A, cs.sparse, cs.iter);  
        hatx = (dctmtx(cs.N))^-1 * recov;
		% cs based csd function
		[Spec, f, alpha] = my_cyclic_spectrum(hatx', cs.N, sig.fs, sig.M,'no-show', 'cs_based_cyclic_spectrum');
		%feature range: f ~[-fs/4 fs/4]; a ~[0, a/4]
		f_mid = round(length(f)/2);
		f_index_range = (f_mid - floor(length(f)/2/3)+1) : (f_mid + floor(length(f)/2/3));
		f_area = f(f_index_range);
		a_index_range = 1 : (ceil(length(alpha) * 0.5));
		a_area = alpha(a_index_range);
		% feature extraction
		y = (Spec(a_index_range,f_index_range));
        energy_fin(kk) = (1/length(cs.N)).*norm(y);% test statistic (normlized)
		%{
		% main component
		Spec_r = reshape(Spec, 1, length(f)* length(alpha));
		sorted_spec = sort(Spec_r, 'descend'); 
		main_spec = sorted_spec(1:10);
		rest_spec = sorted_spec(11:end);
		%test_stat = norm(main_spec) ./ norm(rest_spec);
		%}
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
cs_csd_thresh_est = thresh;
save ./Data/cs_csd_thres.mat cs_csd_thresh_est

