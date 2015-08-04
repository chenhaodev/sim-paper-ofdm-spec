clc; clear; close all

addpath('./Util/')
addpath('./Data/')


% Header 

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
n=(randn(1,sig.N)); 
%sig.x= sig.x + n; 
fs = sig.fs;
N = sig.N;
x = sig.x;
opt1 = 'show';

M = 1;

d_tau = fs/N; % time resolution
tau = 0:d_tau:fs-d_tau; % delay resolution
tau_len = length(tau); 

t_len = floor(N/M-1)+1; 
t_n = -(fs/2-d_tau*floor(M/2)) + d_tau*M*(0:t_len-1); % freq sample location

S = zeros(tau_len, t_len); 
i = 1; 

X = fftshift(fft(x(1:N))); % signal fft
X = X';
x = x';

% Equivalent cyclic spectrum
% Rx = x*x' <=> S in cyclic_spectrum.m 
Dct = dctmtx(N);
Dft = dftmtx(N);
D = Dct; % choose DCT matrix
Rx = x*x';
Sx = D*Rx*D;
%figure; mesh(abs(Sx))

% Vectorize xcorr and cyclic spectrum
matrix_load = 'yes';
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum
Rx_r = reshape(Rx, 1, N*N); %reshape the xcorr
if strcmpi(matrix_load,'yes') % default signal
	load ./Data/matrix_load.mat % load B, H, W_r, H_inv, A 
else
B = eye(N^2);
H = kron((eye(64))',D)*B;
W_r = kron((inv(D))', eye(N));
H_inv = ((H'*H)^(-1))*H';
end
%figure; plot(W_r*Sx_r', 'o'); hold on; plot(H*Rx_r', '*'); %Now H*Rx_r' = W_r*Sx_r';

% Compressed sampling the signal
cs.sparse = 16;
cs.ratio = 8;
cs.iter = 32;
cs.N = N;
cs.M = round(cs.N/cs.ratio); % num of sensing points
Phi = randn(cs.M,cs.N); % sensing (random matrix)
y = Phi*x;

% Link compressed covariance and vectorized cyclic spectrum 
% 1. Phi*Rx*Phi' = Phi*x*x'*Phi' = y*y' = Cy;
% 2. H*Rx_r' = W_r*Sx_r';
% 3. Ry_r = kron(Phi,Phi)*Rx_r'; 
% so Ry_r = kron(Phi,Phi)*H_inv*W_r*Rx_r', where Rx_r is sparse;
Cy = y*y'; %covariance matrix (via compressed data)
Ry_r = reshape(Cy, 1, cs.M*cs.M);
A = kron(Phi,Phi)*H_inv*W_r;
b = Ry_r';

cvx_begin
    variable hatX(N^2);
    minimize(norm(hatX,1));
    A*hatX == b;
cvx_end
hat_m = (vec2mat(hatX, N, N))';
figure; mesh(abs(hat_m));
