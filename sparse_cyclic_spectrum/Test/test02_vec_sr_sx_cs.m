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
%{
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
Dct = dctmtx(N);
Sx = Dct * S * Dct;
%}

% close equivalent cyclic spectrum
Dct = dctmtx(N);
Dft = dftmtx(N);
D = Dct;
Rx = x*x';
Sx = D*Rx*D;
%figure; mesh(abs(Sx))

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

% Dct matrix
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum
Rx_r = reshape(Rx, 1, N*N); %reshape the xcorr (time)
B = eye(N^2);
H = kron((eye(64))',D)*B;
W_r = kron((inv(D))', eye(N));
t1 = W_r*Sx_r';
t2 = H*Rx_r';
%figure; plot(t1, 'o'); hold on; plot(t2, '*');
% Now H*Rx_r' = W_r*Sx_r';
% So Rx_r' = H_inv * W_r*Sx_r'
H_inv = ((H'*H)^(-1))*H';

% compressed xcorr
cs.sparse = 16;
cs.ratio = 8;
cs.iter = 32;
cs.N = N;
cs.M = round(cs.N/cs.ratio);

% sensing (random matrix)
Phi = randn(cs.M,cs.N);
y = Phi*x;

%{
Sy = zeros(cs.M, cs.M); 
i = 1; 

for alfa = tau(1:cs.M)
    interval_tau = round(alfa/d_tau); % tau
    T_N = floor((cs.M-interval_tau-M)/M)+1; 
    
    % window generate
    t = 1:M*T_N;
    t = reshape(t, M, T_N);

    % spectral correlation
	y1 = y(t);
	y2 = y(t+interval_tau); 
   	Syt = conj(y1).*y2;
   	Sy(i, (1:T_N)) = Syt/cs.M; %move St to central
   	i = i+1;
end
%Sy = fftshift(Sy);
%}

%%% In the Real World
Cy = y*y'; %covariance matrix (W/O sum)
% Rx = x*x';   Phi*Rx*Phi' = Phi*x*x'*Phi' = y*y' = Cy;
% since: Rx => S in cyclic_spectrum.m <= Rx
% H*Rx_r' = W_r*Sx_r';
% Ry_r = kron(Phi,Phi)*Rx_r'; 
Ry_r = reshape(Cy, 1, cs.M*cs.M);
A = kron(Phi,Phi)*H_inv*W_r;
b = Ry_r';
%[hat, ~] = cosamp(b, A, 32, 100);

cvx_begin
    variable hatX(N^2);
    minimize(norm(hatX,1));
    A*hatX == b;
cvx_end
hat_m = (vec2mat(hatX, N, N))';
figure; mesh(abs(hat_m));


%{
SY = fftshift(Sy);
Ry_r = reshape(SY, 1, cs.M*cs.M);

PHI = eye(64)*kron(Phi,Phi)*eye(4096);

A = PHI*H_inv*W_r;
b = Ry_r';
%[hat, ~] = cosamp(b, A, 32, 100);
cvx_begin
    variable hatX(4096);
    minimize(norm(hatX,1));
    A*hatX == b;
cvx_end
hat_m = (vec2mat(hatX, 64, 64))';
figure; mesh(abs(hat_m));
%}
