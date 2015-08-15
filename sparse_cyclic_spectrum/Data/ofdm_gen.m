% version 0.1
% the ofdm small system
% read only

% clear
    clear                   % clear all variables
    close all               % close all figures
    clc                     % clear command window
    randn('seed',1214)      % setting the seed for normal random numbers
    rand('seed',12524)      % setting the seed for uniform random numbers
    
% ofdm parameters
    ofdm.N  =  24;                     % number of subcarriers
    ofdm.B  =  1;                       % number of block in each channel realization
    ofdm.M  =  4;                      %ch Modulation order
    ofdm.T  =  1e-7;                    % OFDM sample time
    ofdm.GI =  4;                      % length of gaurd interval
    ofdm.TS =  ofdm.N*ofdm.T;           % OFDM symbol time (not considering gaurd interval)
    ofdm.TT = (ofdm.N+ofdm.GI)*ofdm.T;  % OFDM symbol time (considering gaurd interval)
    ofdm.PP =  1:4:ofdm.N;              % Pilot locations in the subcarriers
    ofdm.DP =  setxor(1:ofdm.N,ofdm.PP);% Data  locations in the subcarriers
    ofdm.NP =  length(ofdm.PP);         % number subcarriers which carry data
    ofdm.groupnum = 5; 
    ofdm.signal = [];
    ofdm.fc = 1;
    ofdm.fs = 3;
    ofdm.rf = 'no';

% config info
    disp('ofdm signal information:'); disp(ofdm);
    
%% Loop
ofdm_signal = [];
for cnt1 = 1 : ofdm.groupnum
        % Data generation
        data  = randi([0 ofdm.M-1],ofdm.N,ofdm.B);
        % modulation
        if ofdm.M == 4
            dataMod = qammod(data,ofdm.M)/sqrt(2);  
        else
            error('Not defined')
        end
        
        % pilot insertion
        ofdm.Pilot = ones(ofdm.NP,1);% or randsrc(ofdm.NP,ofdm.B,[-1 1]).*exp(-sqrt(-1)*pi*rand(ofdm.NP,ofdm.B));
        dataMod(ofdm.PP,:) = ofdm.Pilot; % inseted into 4-QAM dataMod matrix
        
        % ifft operation
        dataIFFT   = sqrt(ofdm.N)*ifft(dataMod); %ifft, now in time domain
        
        % adding gaurd interval
        dataIFFTGI = [dataIFFT((ofdm.N-ofdm.GI+1):ofdm.N)' dataIFFT']'; %copy len=GI dataIFFT as guard interval, insert in front. 
        
        if strcmpi(ofdm.rf, 'yes')
            % RF-front
            dataIFFTGI = dataIFFTGI.';
            dataIFFTGITx = dataIFFTGI.*exp(sqrt(-1)*2*pi*(ofdm.fc)/ofdm.fs*(0:length(dataIFFTGI)-1));
            % RF-receiver
            dataIFFTGIRx = dataIFFTGITx./exp(sqrt(-1)*2*pi*(ofdm.fc)/ofdm.fs*(0:length(dataIFFTGI)-1));
        else
            dataIFFTGITx = dataIFFTGI.';
            dataIFFTGIRx = dataIFFTGITx;
        end
		% OFDM demodulation begin
        % reshaping the signal
        dataChann = reshape(dataIFFTGIRx.',ofdm.N+ofdm.GI,ofdm.B);
        
        % Guard interval removal
        dataRx     = dataChann((ofdm.GI+1):(ofdm.N+ofdm.GI),:);

        % ofdm demodulation
        dataRxFFT  =1/sqrt(ofdm.N)*fft(dataRx);
        dataRxMod_LSE =  dataRxFFT(ofdm.DP,:); % not impulse response, just point scaling  
        dataRxDeMod_LSE = qamdemod(dataRxMod_LSE,ofdm.M);       
        [~,BER_LSE] = biterr(dataRxDeMod_LSE,data(ofdm.DP,:),ofdm.M);
        
        % signal gen
        ofdm.signal = [ofdm.signal dataIFFTGITx];
end
ofdm.signal =  ofdm.signal.*exp(sqrt(-1)*2*pi*(ofdm.fc));

%% Figure
figure; plot(data(ofdm.DP,:), 'o');
hold on; plot(dataRxDeMod_LSE, '*');
save ofdm.mat ofdm
 