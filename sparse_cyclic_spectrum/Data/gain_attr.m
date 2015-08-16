%attr.m 
%attribute for SNR based AWGN signal
%out: signal gain ,noise gain 
%n= randn(1,sig.N); 
%x = sqrt(snr).*s + n;

snr_dB = 10; % SNR in decibels, NOTE !!!!
snr = 10.^(snr_dB./10); % Linear Value of SNR
gain.sig = sqrt(snr);
gain.noise = 1;

save gain_attr.mat gain snr_dB