% ver 1.0, csp based ofdm feature detector
% rayleigh channel model, qam, rand-demodulator (digital and analog)
% xcorr, constellation

addpath('./Util/')

% clear
    clear                   % clear all variables
    close all               % close all figures
    clc                     % clear command window
    randn('seed',1214)      % setting the seed for normal random numbers
    rand('seed',12524)      % setting the seed for uniform random numbers
    
% ofdm parameters
    ofdm.N  =  128;                     % number of subcarriers
    ofdm.B  =  1;                       % number of block in each channel realization
    ofdm.M  =  16;                      %ch Modulation order
    ofdm.T  =  1e-7;                    % OFDM sample time
    ofdm.GI =  16;                      % length of gaurd interval
    ofdm.TS =  ofdm.N*ofdm.T;           % OFDM symbol time (not considering gaurd interval)
    ofdm.TT = (ofdm.N+ofdm.GI)*ofdm.T;  % OFDM symbol time (considering gaurd interval)
    ofdm.PP =  1:10:ofdm.N;              % Pilot locations in the subcarriers
    ofdm.DP =  setxor(1:ofdm.N,ofdm.PP);% Data  locations in the subcarriers
    ofdm.NP =  length(ofdm.PP);         % number subcarriers which carry data
    
% channel parameters
    chan.L      = 10;                               %ch number of channel taps (tab = 1 means only 1 path to RX via channel)
    chan.fd     = .1;                               % doppler in Hz (fD = f*v1/v2 = 2v/lambda)
    chan.Nt     = 128;                              % Number of columns in the dictionary
    chan.Gain   = (0:1/(chan.Nt):1)*0;              % delay spread profile (power sp... pdp?; 129' 0...0)
    [~,chan.Delay]  = sort([0,rand(1,chan.Nt)]);    % generating random delay for each ta[
    chan.snrdBV = 5:2:30;                           % channel signal to noise ration for sweep
    
% loop parameters
    loop.End1  = 1;                               % number of iterations
    loop.End2  = length(chan.snrdBV);               % length of inner loop
    loop.LSE     = zeros(loop.End1,loop.End2);      % memory allocation for the BER using LSE method

% channel type
    chan.choose = 'rayleigh'; % rayleigh or guassian

% compression 
    cs.valid = 'true'; % true or false
    cs.recov = 'false'; % used for further recovery in case of detected ofdm.
    cs.detect = 'true'; 
    cs.digital = 'false';
    ofdm.demod_valid = 'false'; %detection is enough
    
% config info
    
    disp('ofdm signal information:'); disp(ofdm);
    disp('channel information:'); disp(chan);
    disp('compression valid:'); disp(cs);

%% Loop
for cnt1 = 1 :  loop.End1
    for cnt2 = 1 : loop.End2
        % loop parameters
        chan.snrdB = chan.snrdBV(cnt2);
        % Data generation
        data  = randi([0 ofdm.M-1],ofdm.N,ofdm.B);
        % modulation
        if ofdm.M == 16
            dataMod = qammod(data,ofdm.M)/sqrt(2); %4-QAM 
        else
            error('Not defined')
        end
        
        % pilot insertion
        ofdm.Pilot = ones(ofdm.NP,1);% or randsrc(ofdm.NP,ofdm.B,[-1 1]).*exp(-sqrt(-1)*pi*rand(ofdm.NP,ofdm.B));
        dataMod(ofdm.PP,:) = ofdm.Pilot; % in front of 4-QAM dataMod
        %dataMod = zeros(128,1);
        % ifft operation
        dataIFFT   = sqrt(ofdm.N)*ifft(dataMod); %ifft, now in time domain
        
        % adding gaurd interval
        dataIFFTGI = [dataIFFT((ofdm.N-ofdm.GI+1):ofdm.N)' dataIFFT']'; %copy len=GI dataIFFT as guard interval, insert in front. 
        %dataIFFTGI = zeros(144,1);
        
        % channel and transform (rayleigh and gaussian noise)
        if strcmpi(chan.choose,'rayleigh')
            %chan1 (raylei chan) 
            chan.tau_p = linspace(0,ofdm.GI*ofdm.T - ofdm.GI*ofdm.T./chan.Nt,chan.Nt); % each path's delay
            ch = rayleighchan(ofdm.T,chan.fd,chan.tau_p(chan.Delay(1:chan.L)),chan.Gain(chan.Delay(1:chan.L))); % __matlab/rayleighchan.pdf
            dataChann = awgn(filter(ch,dataIFFTGI),chan.snrdB ); %add noise, filter() is considered as convolution, __matlab/filter.rtf
        else
            %chan2 (guassian chan)
            dataChann = awgn(dataIFFTGI,chan.snrdB); %add noise, filter() is considered as convolution, __matlab/filter.rtf
        end
        
        % compression and detection
        if strcmpi(cs.valid,'true')
            cs.iter = 50;
            cs.sparse = 20;
            cs.ratio = 10;
            cs.N = length(dataChann);
            cs.M = round(cs.N/cs.ratio);
            if strcmpi(cs.digital, 'true') % simulate compressed sensing in digital domain 
				cs.psi = (dftmtx(cs.N))^-1;
            	%cs.phi = randn(cs.M,cs.N);   % guassian matrix
            	cs.phi = pn_gen(cs.M,cs.N);   % bernolli matrix
            	cs.y = cs.phi*dataChann; % compressed observation
			else			% simulate compressed sensing in analog domain
				% signal mixing
				pn = pn_gen(cs.N);
				mixed_sig = dataChann.*pn';  
				% ideal pass filter
				lpf_temp = zeros(1,(cs.N-1)); % bandpass: side-band
				lpf_temp(1) = 1;
				lpf_z = interpft(lpf_temp,cs.N)/cs.N; % impulse response
				% low pass filtering and downsampling
				f_lpf = fft(lpf_z');
				f_sig = fft(mixed_sig);
				filter_out = ifft(f_lpf.*f_sig); %multiplication in freq, = convolution in time
				cs.y = downsample(filter_out, cs.ratio); %defactor = downsample ratio
				%cs.temp = downsample(dataChann, cs.ratio);
				%cs.y = cs.temp;
            end
            % feature in cs-sampling: autocorrelation based constellation 
            if strcmpi(cs.detect, 'true')
                figure; plot(xcorr(cs.y),'o'); title('autocorrelated constellation @analog CS-RX(ratio=50), rayleigh chan, 64QAM');
            end
        end
        
        if strcmpi(cs.recov, 'true')
            cs.A = (cs.phi)*(cs.psi); % sensing matrix
            [xhat, trsh] = cosamp(cs.y,cs.A,cs.sparse,cs.iter);
			dataChann = cs.psi*xhat ;
        end 
        
        if strcmpi(ofdm.demod_valid, 'true')
            % OFDM demodulation begin
            dataChann = reshape(dataChann,ofdm.N+ofdm.GI,ofdm.B); % reshaping the signal
            dataRx = dataChann((ofdm.GI+1):(ofdm.N+ofdm.GI),:); % Guard interval removal
            dataRxFFT  =1/sqrt(ofdm.N)*fft(dataRx); % ofdm demodulation
            %% LSE channel estimate 
            H_LSE = zeros(ofdm.N,ofdm.B);
            for b = 1 : ofdm.B
                H_LSE(:,b) = ofdm.N/ofdm.NP * fft(ifft(dataRxFFT(ofdm.PP,b)./dataMod(ofdm.PP,b)),ofdm.N);
            end
            dataRxMod_LSE =  dataRxFFT(ofdm.DP,:)./H_LSE(ofdm.DP,:); % recover signal via chan point response 
            dataRxDeMod_LSE = qamdemod(dataRxMod_LSE,ofdm.M);  % QAM demodulation     
            [~,BER_LSE] = biterr(dataRxDeMod_LSE,data(ofdm.DP,:),ofdm.M);
            loop.LSE(cnt1,cnt2)    = BER_LSE;
        end
    end
    disp([num2str(round(cnt1/loop.End1*100)),'% has been done'])
end

%% Figure
%f1 = figure(1);
%semilogy(chan.snrdBV,mean(loop.LSE,1),'r.-')
%legend('LSE')
%grid on
