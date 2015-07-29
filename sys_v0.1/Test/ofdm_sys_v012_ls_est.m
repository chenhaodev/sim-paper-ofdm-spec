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
    loop.End1  = 1e2;                               % number of iterations
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

% compression 
    cs.valid = 'false'; % true or false

% config info
