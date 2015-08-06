function [hat_m] = sparse_cyclic_spec(x, N, fs, opt1)
% x: signal (1 * N vector)
% N: samples <= len(x) 
% fs: sample rate
% opt1: 'show' => display picture
% author: chenhaomails@gmail.com
% update: 15/08/05

addpath('./Util/')
addpath('./Data/')

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
matrix_load = 'yes'; %NOTE !!! matrix pre-load for faster execution
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum
Rx_r = reshape(Rx, 1, N*N); %reshape the xcorr
if strcmpi(matrix_load,'yes') % default signal
	load ./Data/matrix_load.mat % load B, H, W_r, H_inv, A, Phi, cs (attr)
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
% 3. Ry_r' = kron(Phi,Phi)*Rx_r'; 
% so Ry_r' = kron(Phi,Phi)*H_inv*W_r*Sx_r', where Rx_r is sparse;
% Problem: sparse hat_m is not concentrated!
% Trick: fix the sensing random matrix.
Cy = y*y'; %covariance matrix (via compressed data)
Ry_r = reshape(Cy, 1, cs.M*cs.M);
A = kron(Phi,Phi)*H_inv*W_r;
b = Ry_r';


cvx_begin quiet
    variable hatX(N^2);
    minimize(norm(hatX,1));
    A*hatX == b;
    %subject to
    %    norm( A * hatX - b, 2 ) <= 0.0001;
cvx_end

hat_m = (vec2mat(hatX, N, N))';
if strcmpi(opt1,'show') 
	figure; mesh(abs(hat_m));
end
