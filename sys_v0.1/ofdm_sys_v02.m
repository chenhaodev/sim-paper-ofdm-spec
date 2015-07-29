% ver 0.12, ofdm mini system
% channel model: rayleigh / guassian 
% channel estimate: ls / sparse estimation (training)
% receiver: using the channel estimation
% signal detection: cs-adc + xcorrelation
% display: bit error rate 
% created by chenhaomails@gmail.com

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
    ofdm.M  =  4;                       % Modulation order (M=4)
    ofdm.T  =  1e-7;                    % OFDM sample time
    ofdm.GI =  16;                      % length of gaurd interval, practically equals 4 * ofdm.M
    ofdm.TS =  ofdm.N*ofdm.T;           % OFDM symbol time (not considering gaurd interval)
    ofdm.TT = (ofdm.N+ofdm.GI)*ofdm.T;  % OFDM symbol time (considering gaurd interval)
    ofdm.PP =  1:10:ofdm.N;              % Pilot locations in the subcarriers
    ofdm.DP =  setxor(1:ofdm.N,ofdm.PP);% Data  locations in the subcarriers
    ofdm.NP =  length(ofdm.PP);         % number subcarriers which carry data
    
% channel parameters
    chan.L      = 3;                                % number of channel taps (tab = 1 means only 1 path to RX via channel)
    chan.fd     = .1;                               % doppler in Hz (fD = f*v1/v2 = 2v/lambda)
    chan.Nt     = 128;                              % Number of columns in the dictionary
    chan.Gain   = (0:1/(chan.Nt):1)*0;              % delay spread profile (power sp... pdp?; 129' 0...0)
    [~,chan.Delay]  = sort([0,rand(1,chan.Nt)]);    % generating random delay for each ta[
    chan.snrdBV = 5:2:30;                           % channel signal to noise ration for sweep
    
% loop parameters
	loop.training = 10;
	loop.running = 1;
    loop.End1  = loop.training + loop.running;      % number of iterations
    loop.End2  = length(chan.snrdBV);               % length of inner loop
    loop.Sparse  = zeros(loop.End1,loop.End2);      % memory allocation for the BER using sparse method
    loop.LSE     = zeros(loop.End1,loop.End2);      % memory allocation for the BER using LSE method

% channel type
    chan.choose = 'rayleigh'; % rayleigh or guassian
	chan.est = 'ls' ; % sparse / least square estimation for channel model
    if strcmpi(chan.choose,'rayleigh')
		chan.tau_p = linspace(0,ofdm.GI*ofdm.T - ofdm.GI*ofdm.T./chan.Nt,chan.Nt); % each path's delay
		chan.Gamma = exp(-sqrt(-1)*2*pi.*repmat(((1:ofdm.N).'),1,chan.Nt)./ofdm.TS.*repmat(chan.tau_p,ofdm.N,1));  %1D-dft dictionary with delay info
		disp('other way to build the dict for channel model, see ./chan_est.m')
	end

% compressed sampling 
    cs.detect = 'true'; 
	cs.digital = 'true'; % cs in digital / analog domain
    cs.iter = 50;
    cs.sparse = ofdm.M * 2; % real + image sparsity in phase
    cs.ratio = 10;
    cs.N = ofdm.N;
    cs.M = round(cs.N/cs.ratio);
	cs.psi = sqrt(ofdm.N)* inv(dftmtx(cs.N)); % corres to the ofdm's ifft operation
	dict_4qam.phase = [1 2 3 4]/4;
	dict_4qam.offset = 1/8;
	dict_4qam.D = [exp(2*pi*sqrt(-1)*0) exp(2*pi*sqrt(-1)*dict_4qam.phase + 2*pi*sqrt(-1)*dict_4qam.offset)]; %
    
% config info
    disp('ofdm signal information:'); disp(ofdm);
    disp('channel information:'); disp(chan);
    disp('compression valid:'); disp(cs);
	show.chan_est_ber = 'valid'; 

% Loop ; Traning for channel model
for cnt1 = 1 : loop.End1
    for cnt2 = 1 : loop.End2
        % loop parameters
        chan.snrdB = chan.snrdBV(cnt2);
        % Data generation
        data  = randi([0 ofdm.M-1],ofdm.N,ofdm.B);
        % modulation
        if ofdm.M == 4
            dataMod = qammod(data,ofdm.M)/sqrt(2); %4-QAM 
        else
            error('Not defined')
        end
        % pilot insertion
        ofdm.Pilot = ones(ofdm.NP,1);% or randsrc(ofdm.NP,ofdm.B,[-1 1]).*exp(-sqrt(-1)*pi*rand(ofdm.NP,ofdm.B));
        dataMod(ofdm.PP,:) = ofdm.Pilot; % in front of 4-QAM dataMod
        
        % ifft operation
        dataIFFT   = sqrt(ofdm.N)*ifft(dataMod); %ifft, now in time domain
        
        % adding gaurd interval
        dataIFFTGI = [dataIFFT((ofdm.N-ofdm.GI+1):ofdm.N,:);dataIFFT;]; %copy len=GI dataIFFT as guard interval, insert in front. 
        
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
        
		% OFDM demodulation begin
        % reshaping the signal
        dataChann = reshape(dataChann,ofdm.N+ofdm.GI,ofdm.B);
        
        % Guard interval removal (in practise, windowing)
        dataRx     = dataChann((ofdm.GI+1):(ofdm.N+ofdm.GI),:);

		% low rate detection (without demodu), using compressed sampling
		if strcmpi(cs.detect,'true')
			if (cnt1 > loop.training) 
				pn = pn_gen(cs.N);
				mixed_sig = dataRx.*pn';  % signal mixing
				lpf_temp = zeros(1,(cs.N-1));
				lpf_temp(1) = 1;
				lpf_z = interpft(lpf_temp,cs.N)/cs.N; % ideal pass filter, via interpolation
				f_lpf = fft(lpf_z');
				f_sig = fft(mixed_sig);
				filter_out = ifft(f_lpf.*f_sig); %multiplication in freq, = convolution in time
				cs.y = downsample(filter_out, cs.ratio); % low pass filtering and downsampling
                f1 = figure(1);
				plot(xcorr(cs.y),'o'); title('autocorrelated constellation in cs based rx (ratio=10)'); hold on;
			end
        end

        % ofdm demodulation
        dataRxFFT  =1/sqrt(ofdm.N)*fft(dataRx);
        
		% training step for channel estimation
		if (cnt1 <= loop.training) 
			if strcmpi(chan.est,'sparse') % sparse channel estimation
				H_Sparse = zeros(ofdm.N,ofdm.B);
        		lambda1 = ofdm.NP*10^(-chan.snrdB/10)/sum(abs(ch.pathGains)); %ch.pathGain returned by reyleighchan(),change dB to times,normlize
        		for b = 1 : ofdm.B
        		    y = dataRxFFT(ofdm.PP,b);
        		    A = chan.Gamma(ofdm.PP,:).*repmat(ofdm.Pilot(:,b),1,chan.Nt);
        		    cvx_begin quiet
        		        variable x(chan.Nt) complex
        		            % sparse mini formula (A is built from dictionary, y is received data and x is the chan coeff at pilot locations)
        		            minimize( quad_form(y-A*x,eye(ofdm.NP))+lambda1*norm(x,1) )
        		    cvx_end
        		    % building channel at all location (simply from the dictionary)
        		    H_Sparse(:,b) = chan.Gamma*x;            
        		end
        		dataRxMod_Sparse =  dataRxFFT(ofdm.DP,:)./H_Sparse(ofdm.DP,:); % use estimated chan (by pilots) to derive dataMod         
        		dataRxDeMod_Sparse = qamdemod(dataRxMod_Sparse,ofdm.M);  %  demodu 
        		[~,BER_Sparse] = biterr(dataRxDeMod_Sparse,data(ofdm.DP,:),ofdm.M); % bit error, compare to generated (randi), __matlab/qammod.pdf
				H_est = H_Sparse;
			else % LSE channel estimation        
        		H_LSE = zeros(ofdm.N,ofdm.B);
        		for b = 1 : ofdm.B
        		    H_LSE(:,b) = ofdm.N/ofdm.NP * fft(ifft(dataRxFFT(ofdm.PP,b)./dataMod(ofdm.PP,b)),ofdm.N);
        		end
        		dataRxMod_LSE =  dataRxFFT(ofdm.DP,:)./H_LSE(ofdm.DP,:); % recover signal via chan point response 
        		dataRxDeMod_LSE = qamdemod(dataRxMod_LSE,ofdm.M);  % QAM demodulation     
        		[~,BER_LSE] = biterr(dataRxDeMod_LSE,data(ofdm.DP,:),ofdm.M);
				H_est = H_LSE;
			end
			% saving the output
        	if strcmpi(chan.est,'sparse')
				loop.BER(cnt1,cnt2) = BER_Sparse;
			else
				loop.BER(cnt1,cnt2)    = BER_LSE;
			end
		end 
    end
	if (cnt1 <= loop.training) 
		disp([num2str(round(cnt1/loop.training*100)),'% training has been done'])
	else
		disp([num2str(round((cnt1-loop.training)/loop.running*100)),'% execution has been done'])
    end
    
    % show channel model and biterr
    if (cnt1 == loop.training) 
        if strcmpi(show.chan_est_ber,'valid')
            hold off;
            f2 = figure(2); subplot(2,1,1);
            semilogy(chan.snrdBV,mean(loop.BER,1),'r.-')
            legend(chan.est)
            grid on
            subplot(2,1,2); 
            plot(real(H_est)); title('estimated channel')
        end
    end
   
end

% display
disp('the data input: (e.g. last iter)'); data(ofdm.DP,:)'
disp('the data output: (e.g. last iter)'); dataRxDeMod_LSE'
