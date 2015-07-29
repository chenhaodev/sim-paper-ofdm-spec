% In this mfile it is shown how we can use the sparsity concept for channel
% estimation.

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
    ofdm.GI =  16;                      % length of gaurd interval
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
    loop.End1  = 1e2;                               % number of iterations
    loop.End2  = length(chan.snrdBV);               % length of inner loop
    loop.Sparse  = zeros(loop.End1,loop.End2);      % memory allocation for the BER using sparse method
    loop.LSE     = zeros(loop.End1,loop.End2);      % memory allocation for the BER using LSE method
       
% building dictionary (Gamma, please check different papers to learn how to build the dictionary)
% http://www.mathworks.com/matlabcentral/fileexchange/13127-ofdm-lse-channel-estimation )
    %fft_dict1 chan.Gamma = exp(2*pi*sqrt(-1)/ofdm.N .* meshgrid(0:ofdm.N-1,0:ofdm.N-1).* repmat((0:ofdm.N-1)',[1,ofdm.N])); % fft matrix
    %fft_dict2 chan.Gamma = dftmtx(128);
    %fft_dict3 (http://www.mathworks.com/matlabcentral/fileexchange/39457-ofdm-sparse-channel-estimation), sparse est can be better
    chan.tau_p = linspace(0,ofdm.GI*ofdm.T - ofdm.GI*ofdm.T./chan.Nt,chan.Nt); % each path's delay
    chan.Gamma = exp(-sqrt(-1)*2*pi.*repmat(((1:ofdm.N).'),1,chan.Nt)./ofdm.TS.*repmat(chan.tau_p,ofdm.N,1));  %1D-dft dictionary with delay info
	%fft_dict4 how about overcomplete dict? dft+delay?
    
    %% Loop
for cnt1 = 1 :  loop.End1
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
        ch = rayleighchan(ofdm.T,chan.fd,chan.tau_p(chan.Delay(1:chan.L)),chan.Gain(chan.Delay(1:chan.L))); %create channel, __matlab/rayleighchan.pdf
        dataChann = awgn(filter(ch,dataIFFTGI(:)),chan.snrdB ); %add noise, filter() is considered as convolution, __matlab/filter.rtf
        
		% OFDM demodulation begin
        % reshaping the signal
        dataChann = reshape(dataChann,ofdm.N+ofdm.GI,ofdm.B);
        
        % Guard interval removal
        dataRx     = dataChann((ofdm.GI+1):(ofdm.N+ofdm.GI),:);

        % ofdm demodulation
        dataRxFFT  =1/sqrt(ofdm.N)*fft(dataRx);
        
        %% Saprse Channel estimation
        H_Sparse = zeros(ofdm.N,ofdm.B);
        lambda1 = ofdm.NP*10^(-chan.snrdB/10)/sum(abs(ch.pathGains)); %ch.pathGain is the returned by reyleighchan(), change dB to times, normlize
        for b = 1 : ofdm.B
            y = dataRxFFT(ofdm.PP,b);
            A = chan.Gamma(ofdm.PP,:).*repmat(ofdm.Pilot(:,b),1,chan.Nt);
            cvx_begin quiet
                variable x(chan.Nt) complex
                    % sparse minimization formula (A is built from dictionary, y is received data and x is the channel coeff at pilot locations)
                    minimize( quad_form(y-A*x,eye(ofdm.NP))+lambda1*norm(x,1) )
            cvx_end
            % building channel at all location (simply from the dictionary)
            H_Sparse(:,b) = chan.Gamma*x;            
        end
        dataRxMod_Sparse =  dataRxFFT(ofdm.DP,:)./H_Sparse(ofdm.DP,:); % use estimated chan (by pilots) to derive dataMod         
        dataRxDeMod_Sparse = qamdemod(dataRxMod_Sparse,ofdm.M);  %  demodu 
        [~,BER_Sparse] = biterr(dataRxDeMod_Sparse,data(ofdm.DP,:),ofdm.M); % bit error, compare to generated (randi), __matlab/qammod.pdf

        %% LSE Channel estimation        
        H_LSE = zeros(ofdm.N,ofdm.B);
        for b = 1 : ofdm.B
             H_LSE(:,b) = ofdm.N/ofdm.NP * fft(ifft(dataRxFFT(ofdm.PP,b)./dataMod(ofdm.PP,b)),ofdm.N);
        end
        
        dataRxMod_LSE =  dataRxFFT(ofdm.DP,:)./H_LSE(ofdm.DP,:);             
        dataRxDeMod_LSE = qamdemod(dataRxMod_LSE,ofdm.M);       
        [~,BER_LSE] = biterr(dataRxDeMod_LSE,data(ofdm.DP,:),ofdm.M);
        
        % saving the output
        loop.Sparse(cnt1,cnt2) = BER_Sparse;
        loop.LSE(cnt1,cnt2)    = BER_LSE;
    end
    disp([num2str(round(cnt1/loop.End1*100)),'% has been done'])
end

%% Figure
f1 = figure(1);
semilogy(chan.snrdBV,mean(loop.Sparse,1),'b-o')
hold on 
semilogy(chan.snrdBV,mean(loop.LSE,1),'r.-')
hold off
legend('Sparse','LSE')
grid on
