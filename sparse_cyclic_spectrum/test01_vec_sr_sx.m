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

%% cyclic spec generate

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

% Loop
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
D = dftmtx(N);
W = D*S*D;
Z = fftshift(W);
Spec_f = abs(Z);

% equivalent
%{
for ii = 1:N
    Y(ii,:) = fftshift(fft(S(ii, :)));
end
for jj = 1:N
	Z(:,jj) = fftshift(fft(Y(:, jj)));
end
%}

%% test begin
Sx_r = reshape(Z, 1, 64*64); %reshape the cyclic spectrum
Rx_r = reshape(S, 1, 64*64); %reshape the xcorr (time)
B = eye(4096);
H = kron((eye(64))',D./64)*B;
H_inv = (H)^(-1);
save H_inv_64.mat H_inv
W_r = kron((dftmtx(64))', eye(64));
t1 = H_inv*W_r*Rx_r';
plot(Sx_r, 'o'); hold on; plot(t1, '*'); 
t1_m = (vec2mat(t1, 64, 64))';
figure; mesh(abs(t1_m));


