%This program generate BPSK-OFDM signal.
%No channel and receivers.
%PSD of generated signal is shown.
%Author: chenhaomails@gmail.com

clear all
nFFTSize =64 ;
% for each symbol bits a1 to a52 are assigned to subcarrier 
subcarrierIndex = [-26:-1 1:26]; % index [-26 to -1 1 to 26] 
nBit = 2500; 
ip = rand(1,nBit) > 0.5; % generating 1's and 0's
nBitPerSymbol = 52;

nSymbol = ceil(nBit/nBitPerSymbol);

% BPSK modulation
% bit0 --> -1
% bit1 --> +1
ipMod = 2*ip - 1; %2500    {0 or 1]
ipMod = [ipMod zeros(1,nBitPerSymbol*nSymbol-nBit)]; % 1 x 2548
ipMod = reshape(ipMod,nSymbol,nBitPerSymbol);  % 49 x 52 
figure(1)
plot(ipMod)

st = []; % empty vector

for ii = 1:nSymbol
	inputiFFT = zeros(1,nFFTSize);
	% assigning bits a1 to a52 to subcarriers [-26 to -1, 1 to 26]
	inputiFFT(subcarrierIndex+nFFTSize/2+1) = ipMod(ii,:); 
	%  shift subcarriers at indices [-26 to -1] to fft input indices [38 to 63]
	inputiFFT = fftshift(inputiFFT); 
	% ofdm modulation opera
	outputiFFT = ifft(inputiFFT,nFFTSize);
	% adding cyclic prefix of 16 samples, time domain 
	outputiFFT_with_CP = [outputiFFT(49:64) outputiFFT];
	%transmit 1 ofdm-symbol (with CP) per iter; st is the total trans.
	st = [st outputiFFT_with_CP]; 
end

fsMHz = 20;
[Pxx,W] = pwelch(st,[],[],4096,20); %power density: pwelch(signal,WINDOW,NOVERLAP,NFFT,Fs) 
figure(2)
plot([0:4095]*fsMHz/4095,10*log10(fftshift(Pxx)));
xlabel('frequency, MHz')
ylabel('power spectral density')
title('Transmit spectrum OFDM (based on 802.11a)');
