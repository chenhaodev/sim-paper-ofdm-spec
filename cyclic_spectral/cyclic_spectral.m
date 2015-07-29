% This example illustrates the detection for cyclostationary signal 
% band ~ [.15 .25] 
% cyclic frequency alpha = .00125.
% author: chenhaomails@gmail.com
% 2015.7

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

sig.type = 'default'; % 'default', 'rand', 'ofdm'
sig.snr_dB = 5; % SNR in decibels
sig.cyclic = 16; 

if strcmpi(sig.type,'default') % default signal
	load signal.mat             
elseif strcmpi(sig.type,'rand')% generate cyclic-random signals
	snr = 10.^(sig.snr_dB./10); % Linear Value of SNR
	l = 16;
	L = l * sig.cyclic;
	n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex noise, real + image
	s_temp = sqrt(snr).*(randn(1,l)+sqrt(-1)*randn(1,l))./(sqrt(2));
	s = repmat(s_temp, 1, sig.cyclic); % copy
	x = s + n;
elseif strcmpi(sig.type,'ofdm')% generate cyclic-random signals
	load signal_ofdm.mat 
	[s1, s2, s3] = size(signal_ofdm);
	l = s1;
	L = l * sig.cyclic;
	n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex noise, real + image
	s_temp = reshape(signal_ofdm, 1, s1*s2*s3);
	s = s_temp(1:L);
	x = s + n;
else
	error('signal type not exist!!');
end

cps.L = length(x);		% signal length
cps.Nw = 256;			% window length
cps.Nv = fix(2/3*cps.Nw);	% block overlap
nfft = 2*cps.Nw;		% fft length
cps.da = 1/cps.L;       % cyclic frequency resolution
cps.a1 = 51;            % first cyclic freq. bin to scan
cps.a2 = 200;           % last cyclic freq. bin to scan

% plot the attribution
cps
sig

% Loop over cyclic frequencies
S = zeros(nfft,cps.a2-cps.a1+1);
alpha_set = [cps.a1: cps.da: cps.a2];

for k = cps.a1:cps.a2;
	Sx = cyclic_spectrum(x,x,k/cps.L,nfft,cps.Nv,cps.Nw,'sym'); % Cyclic Power Spectrum vector, CPS.S
	S(:,k-cps.a1+1) = Sx.S; 
end

% Plot results
Fs = 1;			% sampling frequency in Hz
alpha = Fs*(cps.a1:cps.a2)*cps.da;
Coh.f = Sx.f;
f = Fs*Coh.f(1:nfft/2);

figure
subplot(3,1,1);
imagesc(alpha,f,abs(S(1:nfft/2,:)).^2),
colormap(jet),colorbar,axis xy,title('Cyclic Spectral Density'),
xlabel('cyclic frequency \alpha [Hz]'),ylabel('spectral frequency f [Hz]') % observe the significant points, xaxis = 0.00125

subplot(3,1,2);
plot(alpha,abs(S(1:nfft/2,:)).^2);
xlabel('cyclic frequency \alpha [Hz]');

subplot(3,1,3);
plot(f,abs(S(1:nfft/2,:)).^2);
xlabel('spectral frequency f [Hz]');
