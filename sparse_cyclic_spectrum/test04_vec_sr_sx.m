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

X = fftshift(fft(x(1:N))); % signal fft
X = X';
x = x';

% Equivalent cyclic spectrum
% Rx = x*x' <=> S in cyclic_spectrum.m 
Dft = dctmtx(N);%dftmtx(N);
D = fftshift(Dft); % choose DFT matrix
Rx = conj(x)*x';
Sx = conj(D)*Rx*D';

%figure; mesh(abs(Sx))

% Vectorize xcorr and cyclic spectrum
matrix_load = 'yes';
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum
Rx_r = reshape(Rx, 1, N*N); %reshape the xcorr
if strcmpi(matrix_load,'yes') % default signal
	load ./Data/matrix_load_64_8_dct.mat % load D, B, H, W_r, H_inv, cs (attr)
else
%{
tempA = reshape(triu(Rx), 1, N*N);
inxA = find(tempA ~= 0);
Rx_tr = Rx_r(inxA);
B = Rx_r' * Rx_tr ./( norm(Rx_tr)* norm(Rx_tr));
Pn = B;
H = kron((eye(N))',D)*B;
W_r = kron((inv(D))', eye(N));
H_inv = ((H'*H)^(-1))*H';
%}
B = eye(N^2);
H = kron((eye(N))',D)*B;
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
%if strcmpi(matrix_load,'yes') % default signal
%    y = Phi*x;
%else
%Phi = randn(cs.M,cs.N); % sensing (random matrix)
Phi = pn_gen(cs.M,cs.N);
%end
y = Phi*x;


% Link compressed covariance and vectorized cyclic spectrum 
% 1 Phi*conj(x)*x'*Phi' = Phi*Rx*Phi' =conj(y)*y' = Cy; ..ok
% 2. H*Rx_r' = W_r*Sx_r'; ..ok
% 3. Cy_r' = kron(Phi,Phi)*Rx_r'; .. ok
% so Cy_r' = kron(Phi,Phi)*H_inv*W_r*Rx_r', where Rx_r is sparse;
% Problem: sparse hat_m is not concentrated!
%Cy = conj(y)*y'; %covariance matrix (via compressed data)

%{
Cy = y*y'; %covariance matrix (via compressed data)
Cy_r = reshape(Cy, 1, cs.M*cs.M);
tempB = reshape(triu(Cy), 1, cs.M*cs.M);
inxB = find(tempB ~= 0);
Cy_tr = Cy_r(inxB);
Qm = Cy_tr' * Cy_r ./( norm(Cy_r)* norm(Cy_r));
A = Qm*kron(Phi,Phi)*Pn*H_inv*W_r;
b = Cy_tr';
%}
Cy = conj(y)*y'; %covariance matrix (via compressed data)
Cy_r = reshape(Cy, 1, cs.M*cs.M);
K = kron(D', eye(N));
A = kron(Phi,Phi)*H_inv*W_r*K;
b = Cy_r';

cvx_begin
    variable hatX(N^2);
    minimize(norm(hatX,1));
    A*hatX == b;
cvx_end
hat_m = (vec2mat(hatX, N, N))';
figure; mesh(abs(hat_m));
figure; mesh(abs(hat_m(1:N/2, 1:N/2)));

% Extract feature energy from recov spectrum
%[out] = feature_extract(abs(hat_m), 1:N, 0.2, 1:N, 0.2);
%norm(out)
