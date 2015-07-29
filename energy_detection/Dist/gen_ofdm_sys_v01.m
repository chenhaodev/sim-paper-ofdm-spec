% version 0.1
% the ofdm small system
% generate ofdm signal
% 

% clear
    clear                   % clear all variables
    close all               % close all figures
    clc                     % clear command window
	addpath('../Util')
    randn('seed',1214)      % setting the seed for normal random numbers
    rand('seed',12524)      % setting the seed for uniform random numbers
    
% ofdm parameters
    ofdm.N  =  128;                     % number of subcarriers
    ofdm.B  =  1;                       % number of block in each channel realization
    ofdm.M  =  4;                      %ch Modulation order
    ofdm.T  =  1e-7;                    % OFDM sample time
    ofdm.GI =  16;                      % length of gaurd interval
    ofdm.TS =  ofdm.N*ofdm.T;           % OFDM symbol time (not considering gaurd interval)
    ofdm.TT = (ofdm.N+ofdm.GI)*ofdm.T;  % OFDM symbol time (considering gaurd interval)
    ofdm.PP =  1:10:ofdm.N;              % Pilot locations in the subcarriers
    ofdm.DP =  setxor(1:ofdm.N,ofdm.PP);% Data  locations in the subcarriers
    ofdm.NP =  length(ofdm.PP);         % number subcarriers which carry data (valid subcarriers)
    
% channel parameters
    chan.L      = 3;                               %ch number of channel taps (tab = 1 means only 1 path to RX via channel)
    chan.fd     = .1;                               % doppler in Hz (fD = f*v1/v2 = 2v/lambda)
    chan.Nt     = 128;                              % Number of columns in the dictionary
    chan.Gain   = (0:1/(chan.Nt):1)*0;              % delay spread profile (power sp... pdp?; 129' 0...0)
    [~,chan.Delay]  = sort([0,rand(1,chan.Nt)]);    % generating random delay for each ta[
%    chan.snrdBV = 5:2:30;                           % channel signal to noise ration for sweep
    chan.snrdBV = 0;                           % channel signal to noise ration for sweep
    
% loop parameters
    loop.End1  = 100;                               % number of iterations
    loop.End2  = length(chan.snrdBV);               % length of inner loop
    loop.LSE     = zeros(loop.End1,loop.End2);      % memory allocation for the BER using LSE method
       
% building dictionary (Gamma, please check different papers to learn how to build the dictionary)
    chan.tau_p = linspace(0,ofdm.GI*ofdm.T - ofdm.GI*ofdm.T./chan.Nt,chan.Nt); % each path's delay

% config info
    disp('ofdm signal information:'); disp(ofdm);
    disp('channel information:'); disp(chan);
    disp('loop parameters:'); disp(loop);
    
%% Loop
for cnt1 = 1 :  loop.End1
    for cnt2 = 1 : loop.End2
        % loop parameters
        chan.snrdB = chan.snrdBV(cnt2);
        % Data generation
        data  = randi([0 ofdm.M-1],ofdm.N,ofdm.B);
        %mask  = randi([0 1],1,ofdm.B);
        % modulation
        if ofdm.M == 4
            dataMod = qammod(data,ofdm.M)/sqrt(2);  
        else
            error('Not defined')
        end
        
        % pilot insertion
        ofdm.Pilot = ones(ofdm.NP,1);% or randsrc(ofdm.NP,ofdm.B,[-1 1]).*exp(-sqrt(-1)*pi*rand(ofdm.NP,ofdm.B));
        dataMod(ofdm.PP,:) = ofdm.Pilot; % in front of 4-QAM dataMod
        
        % ifft operation
        dataIFFT   = sqrt(ofdm.N)*ifft(dataMod); %ifft, now in time domain
        
        % adding gaurd interval
        dataIFFTGI = [dataIFFT((ofdm.N-ofdm.GI+1):ofdm.N)' dataIFFT']'; %copy len=GI dataIFFT as guard interval, insert in front. 
		%if mask == 0
		%	dataIFFTGI = zeros(length(dataIFFTGI),1);
		%end

        % channel and transform (rayleigh and gaussian noise)
        ch = rayleighchan(ofdm.T,chan.fd,chan.tau_p(chan.Delay(1:chan.L)),chan.Gain(chan.Delay(1:chan.L))); %create channel, __matlab/rayleighchan.pdf
        dataChann = awgn(filter(ch,dataIFFTGI),chan.snrdB ); %add noise, filter() is considered as convolution, __matlab/filter.rtf
		
		% generated ofdm siganl
		signal_ofdm(:, cnt2, cnt1) = dataChann;
		%signal_mask(:, cnt2, cnt1) = mask;
		%signal_data(:, cnt2, cnt1) = data;

    %{    
        
		% OFDM demodulation begin
        % reshaping the signal
        dataChann = reshape(dataChann,ofdm.N+ofdm.GI,ofdm.B);
        
        % Guard interval removal
        dataRx     = dataChann((ofdm.GI+1):(ofdm.N+ofdm.GI),:);

        % ofdm demodulation
        dataRxFFT  =1/sqrt(ofdm.N)*fft(dataRx);
        %% LSE        
        H_LSE = zeros(ofdm.N,ofdm.B);
        for b = 1 : ofdm.B
             H_LSE(:,b) = ofdm.N/ofdm.NP * fft(ifft(dataRxFFT(ofdm.PP,b)./dataMod(ofdm.PP,b)),ofdm.N); % point divid -> point multiply
        end

        dataRxMod_LSE =  dataRxFFT(ofdm.DP,:)./H_LSE(ofdm.DP,:); % not impulse response, just point scaling  
        dataRxDeMod_LSE = qamdemod(dataRxMod_LSE,ofdm.M);       
        [~,BER_LSE] = biterr(dataRxDeMod_LSE,data(ofdm.DP,:),ofdm.M);
        loop.LSE(cnt1,cnt2)    = BER_LSE;
	%}
    end
    disp([num2str(round(cnt1/loop.End1*100)),'% has been done'])
end
save Data/sig_ofdm/signal_ofdm.mat signal_ofdm
%save Data/sig_ofdm/signal_mask.mat signal_mask
%save Data/sig_ofdm/signal_data.mat signal_data

%% Figure
%f1 = figure(1);
%semilogy(chan.snrdBV,mean(loop.LSE,1),'r.-')
%legend('LSE')
%grid on
