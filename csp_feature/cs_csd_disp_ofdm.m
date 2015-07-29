% This code demo the cyclic spectrum density of cyclostationary signal (e.g. ofdm) 
% author: chenhaomails@gmail.com
% 2015.7

clc; clear; close all

addpath('../Util/')
addpath('../Data/')

sig.type = 'ofdm_auto'; % 'default', 'rand', 'ofdm'
sig.snr_dB = 10; % SNR in decibels
sig.cyclic = 16; 
sig.fs = 4; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 10;
snr = 10.^(sig.snr_dB./10); % Linear Value of SNR


if strcmpi(sig.type,'default') % default signal
	load signal.mat            
    x = x(1:2048);
elseif strcmpi(sig.type,'rand')% generate cyclic-random signals
	l = 128;
	L = l * sig.cyclic;
	n=(randn(1,L)+sqrt(-1)*randn(1,L))./(sqrt(2)); % complex noise, real + image
	s_temp = sqrt(snr).*(randn(1,l)+sqrt(-1)*randn(1,l))./(sqrt(2));
	s = repmat(s_temp, 1, sig.cyclic); % copy
	x = s + n;
elseif strcmpi(sig.type,'qam')% generate qam signals
    load ofdm_attr.mat
	data  = randi([0 ofdm.M-1],1,ofdm.N * sig.cyclic);
    x = qammod(data,ofdm.M)/sqrt(2); 
    x = interp(x, (sig.fs)/(sig.fc)); %upsample
    %x = x.*exp(sqrt(-1)*2*pi*(sig.fc)/(sig.fs)*(0:length(x)-1));
    x = awgn(x, sig.snr_dB);
elseif strcmpi(sig.type,'ofdm')% generate ofdm signals
    load ofdm_attr.mat
    for kk = 1:1:sig.cyclic
        data  = randi([0 ofdm.M-1],ofdm.N,ofdm.B);
        dataMod = qammod(data,ofdm.M)/sqrt(2); 
        ofdm.Pilot = ones(ofdm.NP,1);
        dataMod(ofdm.PP,:) = ofdm.Pilot; 
        dataIFFT   = sqrt(ofdm.N)*ifft(dataMod);
        dataIFFTGI = [dataIFFT((ofdm.N-ofdm.GI+1):ofdm.N)' dataIFFT']';
        x(:,kk) = dataIFFTGI;
    end
    x = reshape(x, 1, size(x,1)*sig.cyclic);
    x = interp(x, (sig.fs)/(sig.fc)); %upsample
    %x = x.*exp(sqrt(-1)*2*pi*(sig.fc)/(sig.fs)*(0:length(x)-1));
    x = awgn(x, sig.snr_dB);
    %x = awgn(zeros(1,length(x)), sig.snr_dB); %noise
    csd_sig = x;
else
    load signal_ofdm.mat
    [s1, s2, s3] = size(signal_ofdm);
    %x_temp = signal_ofdm(:,1,(1:sig.cyclic));
    %x = reshape(x_temp, 1, s1*sig.cyclic);
    x = signal_ofdm(:,1,1)';
    %x =(randn(1,144)+sqrt(-1)*randn(1,144))./(sqrt(2)); 
end

csd_sig = x;

% signal length
cps.L = length(x);		

% compression 
cs.sparse = 32; % real + image sparsity in phase
cs.ratio = 4;
cs.N = cps.L;
cs.M = round(cs.N/cs.ratio);
% cs sampling 
pn = pn_gen(cs.N);
mixed_sig = x.*pn;  % signal mixing
lpf_temp = zeros(1,(cs.N-1));
lpf_temp(1) = 1;
lpf_z = interpft(lpf_temp,cs.N)/cs.N; % ideal pass filter, via interpolation
f_lpf = fft(lpf_z);
f_sig = fft(mixed_sig);
filter_out = ifft(f_lpf.*f_sig); %multiplication in freq, = convolution in time
cs.y = downsample(filter_out, cs.ratio); % low pass filtering and downsampling

[Spec, f, alpha] = cyclic_spectrum_new(cs.y, length(cs.y), sig.fs/sig.fc, sig.M, 'show');