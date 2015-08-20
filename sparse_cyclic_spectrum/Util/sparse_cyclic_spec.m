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

% Vectorize xcorr and cyclic spectrum
matrix_load = 'no'; %NOTE !!! matrix pre-load for faster execution
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum
Rx_r = reshape(Rx, 1, N*N); %reshape the xcorr
if strcmpi(matrix_load,'yes') % default signal
	load ./Data/matrix_load.mat % load B, H, W_r, H_inv, A, Phi, cs (attr)
else
B = eye(N^2);
H = kron((eye(N))',D)*B;
W_r = kron((inv(D))', eye(N));
H_inv = ((H'*H)^(-1))*H';
end

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
cs.ratio = 8;
cs.iter = 32;
cs.N = N;
cs.M = round(cs.N/cs.ratio); % num of sensing points
Phi = randn(cs.M,cs.N); % sensing (random matrix)
y = Phi*x;
Cy = y*y'; %covariance matrix (via compressed data)
Ry_r = reshape(Cy, 1, cs.M*cs.M);
A = kron(Phi,Phi)*H_inv*W_r; % equivalent sensing matrix for CS
b = Ry_r';

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
