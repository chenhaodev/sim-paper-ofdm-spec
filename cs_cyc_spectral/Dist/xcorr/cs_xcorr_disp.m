% This code demo the cyclic spectrum density of cyclostationary signal (e.g. ofdm) 
% author: chenhaomails@gmail.com
% 2015.7

clc; clear; close all

addpath('../../Util/')
addpath('../../Data/')

sig.type = 'fsk'; % 'default', 'rand', 'ofdm'
sig.snr_dB = 10; % SNR in decibels
sig.cyclic = 16; 
sig.fs = 4; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 2;
snr = 10.^(sig.snr_dB./10); % Linear Value of SNR

if strcmpi(sig.type,'default') % default signal
	load signal.mat            
    x = x(1:2048);
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
    x = x.*exp(sqrt(-1)*2*pi*(sig.fc)/(sig.fs)*(0:length(x)-1));
    x = awgn(x, sig.snr_dB); %real
    %x = awgn(zeros(1,length(x)), sig.snr_dB); %noise
else
    load FSK_real.mat
    %load signal_ofdm.mat
    x = fsk_real(1:512);
	%x =(randn(1,512)+sqrt(-1)*randn(1,512))./(sqrt(2));
	x = x./norm(x);
end
csd_sig = x;

% autocorrelation matrix
rxx = acf_mtx(x); 

% cyclic spectrum detection
cps.L = length(x);		
[Spec, f, alpha] = cyclic_spectrum_new(csd_sig, cps.L, sig.fs, sig.M,'no-disp');

% dictionary
spec = [Spec Spec];
D = spec * inv(rxx) ; % spec = D * rxx; then rxx = inv(D)* spec, spec is sparse.

% compression 
cs.sparse = 16; % real + image sparsity in phase
cs.ratio = 16;
cs.iter = 50;
cs.N = length(x);
cs.M = round(cs.N/cs.ratio);
% cs sampling 
Phi = randn(cs.M,cs.N);
y = Phi*x';
rxx_y = acf_mtx(y);
A = Phi*D;
for ii = 1:cs.M
    [recov1,~] = cosamp(rxx_y(:,ii), A, cs.sparse, cs.iter); % recover vec(spec) 
    xhat1(:,ii) = recov1;
end

for jj = 1:cs.M
    [recov2,~] = cosamp(rxx_y(jj,:)', A, cs.sparse, cs.iter); % recover vec(spec) 
    xhat2(jj,:) = recov2;
end

figure; mesh(1:cs.M, alpha, abs(xhat1));
xlabel('1:cs.M'); ylabel('alpha');   

figure; mesh([f f], 1:cs.M, abs(xhat2));
xlabel('f'); ylabel('1:cs.M');   


