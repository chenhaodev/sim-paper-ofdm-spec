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
%sig.x = fsk_real(65:128);
%sig.x=(randn(1,64)); 

sig.x = sig.x ./ norm(sig.x);
sig.N = length(sig.x);
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

x = x';

% Equivalent cyclic spectrum
% Rx = x*x' <=> S in cyclic_spectrum.m 
Dft = dctmtx(N);%dftmtx(N);
%D = fftshift(Dft); % choose DFT matrix
D = Dft;
Rx = x*x';
Sx = D*Rx*D;
% add your new idea bgn
D = Dft;
Rx = x*x';
%generate cyclic spectrum 
%ver1: Sx = D*Rx*D;
%ver2:
Sx_tmp = zeros(N,N);
for kk = 1:N 
	Eye_tmp = zeros(N,N);
	Eye_tmp(kk,kk) = 1; % make groups of diag matrix that all entries are '0' but (kk,kk) = 1
	Sx_tmp = D*Rx*Eye_tmp; 
	Sx_tmp_sum = Sx_tmp_sum + Sx_tmp;
end
Sx = Sx_tmp_sum*D;

% add your new idea end

% Vectorize xcorr and cyclic spectrum
matrix_load = 'no';
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum
Rx_r = reshape(Rx, 1, N*N); %reshape the xcorr
B = eye(N^2);
%H = kron((eye(N))',D)*B;
tmp_mtx = zeros(N^2,N^2);
for kk = 1:N % calculate the matrix that map Rx => Sx_tmp_sum
	Eye_tmp = zeros(N,N);
	Eye_tmp(kk,kk) = 1; 
	tmp_mtx = tmp_mtx + kron(D', Eye_tmp); 
end
H = tmp_mtx* B; %something wrong here, since Eye_tmp(k,k) should correspond to G(k), while here we use D instead.
W_r = kron((inv(D))', eye(N));
H_inv = ((H'*H)^(-1))*H'; %rank(H) = 4096;

% assert: 
t1a = W_r*Sx_r';
t1b = H*Rx_r';
if ( (norm(t1a - t1b) / length(t1a)) < 0.0001)
    disp('test: H*Rx_r == W_r*Sx_r (?) ... yes');
else
    error('test: H*Rx_r == W_r*Sx_r (?) ... no');
end
t2a = H_inv*W_r*Sx_r';
t2b = Rx_r';
if ( (norm(t2a - t2b) / length(t2a)) < 0.0001)
    disp('test: Rx_r == H_inv*W_r*Sx_r (?) ... yes');
else
    error('test: Rx_r == H_inv*W_r*Sx_r (?) ... no');
end

% Compressed sampling the signal
cs.sparse = 16;
cs.ratio = 4;
cs.iter = 32;
cs.N = N;
cs.M = round(cs.N/cs.ratio); % num of sensing points
Phi = pn_gen(cs.M,cs.N);
y = Phi*x;
Cy = y*y'; %covariance matrix (via compressed data)
Cy_r = reshape(Cy, 1, cs.M*cs.M);
A = kron(Phi,Phi)*H_inv*W_r; % equivalent sensing matrix for CS
b = Cy_r';

% assert: 
t3a = A*Sx_r';
t3b = b;
if ( (norm(t3a - t3b) / length(t3a)) < 0.0001)
    disp('test: Cy_r = A*Sx_r (?) ... yes');
else
    error('test: Cy_r = A*Sx_r (?) ... no');
end

% Link compressed covariance and vectorized cyclic spectrum 
% 1 Phi*x*x'*Phi' = Phi*Rx*Phi' =y*y' = Cy; ..ok
% 2. H*Rx_r' = W_r*Sx_r'; ..ok
% 3. Cy_r' = kron(Phi,Phi)*Rx_r'; .. ok
% so Cy_r' = kron(Phi,Phi)*H_inv*W_r*Rx_r', where Rx_r is sparse;
% Problem: sparse hat_m is not concentrated!

cvx_begin
    variable hatX(N^2);
    minimize(norm(hatX,1));
    A*hatX == b;
cvx_end
hat_m = (vec2mat(hatX, N, N))';
figure; mesh(abs(hat_m));

% Extract feature energy from recov spectrum
%[out] = feature_extract(abs(hat_m), 1:N, 0.2, 1:N, 0.2);
%norm(out)
