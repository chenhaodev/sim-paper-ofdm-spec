% calculate the threhsold in signal detection
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
fs = 1;  % normlized sampling rate
fc = 2; %normlized RF-carrier rate
M = 2; % window length
iter = 10;
sig.type = 'fsk'; % 'default', 'rand', 'ofdm'
%sig.snr_dB = 10; % SNR in decibels
sig.fs = 1; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 2;

load FSK_real.mat
ref = fsk_real(1:256);
% signal length
cps.L = length(ref);	

% resolution 
d_alpha = sig.fs/cps.L; % freq resolution
alpha = 0:d_alpha:(sig.fs-d_alpha); % cyclic resolution
a_len = length(alpha); 
f_len = cps.L; 

%% loop %%

for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise
		n=(randn(1,cps.L)+sqrt(-1)*randn(1,cps.L))./(sqrt(2)); % complex, real + image
		% input
        x = sqrt(snr).*ref + n; 
		% cs cyclic spectrum detection
		[Spec, f, alpha] = cyclic_spectrum_new(x, cps.L, sig.fs, sig.M,'no-disp');
		% dictionary
		rxx = acf_mtx(x); 
		spec = [Spec Spec];
		D = spec * inv(rxx) ; % spec = D * rxx; then rxx = inv(D)* spec, spec is sparse.
		% compression 
		cs.sparse = 8; % real + image sparsity in phase
		cs.ratio = 8;
		cs.iter = 10;
		cs.N = length(x);
		cs.M = round(cs.N/cs.ratio);
		% cs sampling 
		Phi = randn(cs.M,cs.N);
		y = Phi*x';
		rxx_y = acf_mtx(y);
		A = Phi*D;
		% cs recov
		for ii = 1:cs.M
		    [recov1,~] = cosamp(rxx_y(:,ii), A, cs.sparse, cs.iter); % recover vec(spec) 
		    xhat1(:,ii) = recov1;
		end
		for jj = 1:cs.M
		    [recov2,~] = cosamp(rxx_y(jj,:)', A, cs.sparse, cs.iter); % recover vec(spec) 
		    xhat2(jj,:) = recov2;
        end
        test_1 = xhat1(round(cs.N * 0.2): round(cs.N * 0.3), :);
 		test_2 = xhat2(: ,  round(cs.N * 0.45): round(cs.N * 0.55));
		test_statistic = norm(test_1) + norm(test_2);
		%energy = abs(y).^2;  % feature-power
		%energy_fin(kk) =(1/L).*sum(energy); % test statistic (normlized)
        energy_fin(kk) = (1/length(cps.L)).*test_statistic;% feature-eig, test statistic (normlized)
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
cs_xcorr_fsk_thresh_est = thresh;
save ./Data/cs_xcorr_fsk_thres.mat cs_xcorr_fsk_thresh_est

