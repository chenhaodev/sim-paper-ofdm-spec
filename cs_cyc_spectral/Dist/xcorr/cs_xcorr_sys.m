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
fs = 1;  % normlized sampling rate
fc = 2; %normlized RF-carrier rate
M = 2; % window length
iter = 20;
sig.type = 'fsk'; % 'default', 'rand', 'ofdm'
%sig.snr_dB = 10; % SNR in decibels
sig.fs = 1; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 2;
load ./Data/cs_xcorr_fsk_thres.mat %cs_xcorr_fsk_thresh_est

load FSK_real.mat
ref = fsk_real(1:256);
% signal length
cps.L = length(ref);	

for m = 1:length(Pf)
    m
    i = 0;
	for kk=1:iter % Number of Monte Carlo Simulations
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
        energy_fin = (1/length(cps.L)).*test_statistic;% feature-eig, test statistic (normlized)
		% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
		thresh(m) = cs_xcorr_fsk_thresh_est(m); %use estimated threshold to detect ofdm
		if(energy_fin >= thresh(m))  
			i = i+1;
		end
    end
    Pd_cs_xcorr_th_est(m) = i/kk; 
end
plot(Pf, Pd_cs_xcorr_th_est, '*r'); hold on; 
save ./Data/cs_xcorr_sys.mat Pd_cs_xcorr_th_est 
load ./Data/ed_ofdm_sys.mat
plot(Pf, Pd_rnd, 'ob'); hold on;
plot(Pf, Pd_ofdm, '.y'); hold on;
plot(Pf, Pd_ofdm_th_est, '.g'); hold on;
legend('ROC of FD for FSK using estimated threshold','ROC of ED for random signal', 'ROC of ED for OFDM', 'ROC of ED for OFDM using estimated threshold');
xlabel('Pf');
ylabel('Pd');


