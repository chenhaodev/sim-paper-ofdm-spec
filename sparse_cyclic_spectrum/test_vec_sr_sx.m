clc; clear; close all

addpath('./Util/')
addpath('./Data/')

sig.type = 'fsk'; % 'fsk'
sig.fs = 1;
sig.M = 1;

if strcmpi(sig.type,'fsk') % default signal
	load fsk.mat             
else
	error('signal type not exist!!');
end

sig.x = fsk_real(1:64);
sig.N = length(sig.x);

fs = sig.fs;
N = sig.N;
x = sig.x;
opt1 = 'show';
%[Spec_cs, Spec, f, a] = cs_cyclic_spectrum(sig.x, sig.N, sig.fs, 'show'); 
%function [Spec_f_cs, Spec_f, f, alpha] = cs_cyclic_spectrum(x, N, fs, opt1)
% x: signal (1 * N vector)
% N: samples <= len(x) 
% fs: sample rate
% author: chenhaomails@gmail.com
% opt1: 'show' => display picture

M = 1;

d_tau = fs/N; % time resolution
tau = 0:d_tau:fs-d_tau; % delay resolution
tau_len = length(tau); 

t_len = floor(N/M-1)+1; 
t_n = -(fs/2-d_tau*floor(M/2)) + d_tau*M*(0:t_len-1); % freq sample location

S = zeros(tau_len, t_len); 
i = 1; 

% signal fft
X = fftshift(fft(x(1:N))); 
X = X';
x = x';

%normal xcorr and cyclic_spec

%% Loop
for alfa = tau

    interval_tau = round(alfa/d_tau); % tau
    T_N = floor((N-interval_tau-M)/M)+1; 
    
    t = 1:M*T_N;
    t = reshape(t, M, T_N);

    % spectral correlation
	x1 = x(t);
	x2 = x(t+interval_tau); 
   	St = conj(x1).*x2;
   	S(i, floor((t_len-T_N)/2)+(1:T_N)) = St/N; %move St to central
   	i = i+1;
end

Spec_t = abs(S);
%for ii = 1:N
%	%Y(ii,:) = simple_fft_tau(S(ii, :), N, tau);
%    Y(ii,:) = fftshift(fft(S(ii, :)));
%end
%for jj = 1:N
%	Z(:,jj) = fftshift(fft(Y(:, jj)));
%end
D = dftmtx(N);
W = D*S*D;
Z = fftshift(W);
Spec_f = abs(Z);


%% test begin
Sf_r = reshape(Z, 1, 64*64);
S_r = reshape(S, 1, 64*64);
B_r = eye(4096);
H_r = kron((eye(64))',D./64)*B_r;
H_inv_r = (H_r)^(-1);
save H_inv_r_64.mat H_inv_r
D_r = kron((dftmtx(64))', eye(64));
t1 = H_inv_r*D_r*S_r';
plot(Sf_r, 'o'); hold on; plot(t1, '*');