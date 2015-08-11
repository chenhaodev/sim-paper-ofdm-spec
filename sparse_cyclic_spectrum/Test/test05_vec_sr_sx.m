% NOTE: There's tricks in generating the Qm, Pn matrix

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

x = x.';

% Equivalent cyclic spectrum
Dft = dctmtx(N);%dftmtx(N);
%D = fftshift(Dft); % choose DFT matrix
D = Dft;

%the following canbe optimised !!!

%R generation, ref[1].eq(7)
R = zeros(N,N);
for nn = 1:N
	for v = 1:(N-nn+1) % last nn-1 items in each row is zero.
        R(nn,v) = x(nn)* x(nn+v-1);
	end
end

%vec{R} generate, ref[1].eq(8)
B = [];
tmp = eye(N*N);
for nn = 1:N
	for v = 1:(N-nn+1) 
		q = (v-1)*N + (nn-1) + 1;  %q = nN + v, ~[0, N-1] => ~[1,N]
		B = [B, tmp(:, q)];
	end
end

%rx generate, ref[1].eq(8)
rx = [];
for nn = 1:N
	for v = 1:(N-nn+1) 
		tmp = R(nn,v);
		rx = [rx, tmp];
	end
end

% assert: 
t0a = B*rx.'; 
t0b = reshape(R, N*N, 1);
if ( (norm(t0a - t0b) / length(t0a)) < 1e-10)
    disp(' vec{R} = B*rx (?) ... yes');
else
    error(' vec{R} = B*rx (?) ... no');
end

%generate cyclic spectrum 
Rxc = zeros(N,N);
for v = 1:N 
	Gv = zeros(N,N);
	Dv = zeros(N,N);
	for nn = 1:N
		for alph = 1:N
			%Gv(nn, alph) = (1/N)*exp(-2*pi*sqrt(-1)*(alph-1)*(nn-1+((v-1)/2))/N); % recall dftmtx(N): each W(j,k)= exp(-2*pi*sqrt(-1)*(j-1)*(k-1)/N); 
			Gv(alph, nn) = (1/N)*exp(-2*pi*sqrt(-1)*(alph-1)*(nn-1+((v-1)/2))/N); % recall dftmtx(N): each W(j,k)= exp(-2*pi*sqrt(-1)*(j-1)*(k-1)/N); 
		end
	end
	Dv(v,v) = 1;  % make groups of diag matrix that all entries are '0' but (v,v) = 1
	Rxc = Rxc + Gv*R*Dv; % ref[1].eq(9)
	% save Gv, Dv
	Gv_save(:,:,v) = Gv;
	Dv_save(:,:,v) = Dv;
end
Sx = Rxc*D; %ref[1].eq(10)

% Vectorize xcorr and cyclic spectrum
matrix_load = 'no';
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum
H_tmp = zeros(N^2,N^2);
for v = 1:N % calculate the matrix that map R => Rxc_sum
	Gv = Gv_save(:,:,v);
	Dv = Dv_save(:,:,v);
	H_tmp = H_tmp + kron(Dv.', Gv);  
end
H = H_tmp* B; 
W_r = kron((inv(D)).', eye(N));
H_inv = ((H'*H)^(-1))*H'; %H*

% assert: 
% pinv(W_r)*H*rx.' = Sx_r.';
t1a = W_r*Sx_r.'; 
t1b = H*rx.';
if ( (norm(t1a - t1b) / length(t1a)) < 1e-10)
    disp('test: H*rx == W_r*Sx_r (?) ... yes');
else
    error('test: H*rx == W_r*Sx_r (?) ... no');
end
t2a = H_inv*W_r*Sx_r.';
t2b = rx.';
if ( (norm(t2a - t2b) / length(t2a)) < 1e-10) % problem, not pass, why ?
    disp('test: rx == H_inv*W_r*Sx_r (?) ... yes');
else
    error('test: rx == H_inv*W_r*Sx_r (?) ... no');
end

% Compressed sampling the signal
cs.sparse = 16;
cs.ratio = 4;
cs.iter = 32;
cs.N = N;
cs.M = round(cs.N/cs.ratio); % num of sensing points
M = cs.M;
Phi = pn_gen(cs.M,cs.N);
y = Phi*x;
Rz = y*y.';

%vec{Rz}
rz = [];
for nn = 1:M
	for v = 1:(M-nn+1) 
		tmp = Rz(nn,v);
		rz = [rz, tmp];
	end
end

%Pn, Qm generate, ref[1].eq(43,44)
%{
Pn = zeros(N*N, N*(N+1)/2);
for v = 1:N
	for nn = 1:(N-v+1)  %eq(43)
		p = (v-1)*N - (v-1)*(v-2)/2 + (nn-1) + 1;
		q1 = (v-1+nn-1)*N + (nn-1) + 1;  %eq(42)
		q2 = (nn-1)*N + (nn-1) + (v-1) + 1;  
		Pn(q1, p) = 1; Pn(q2, p) = 1;
	end
end
%}

Pn = zeros(N*N, N*(N+1)/2);
for v = 0:N-1
	for nn = 0:(N-v-1)  %eq(43)
		p = v*N - (v-1)*v/2 + nn ;
		q1 = (v+nn)*N + nn;  %eq(42)
		q2 = nn*N + nn + v;  
		Pn(q1+1, p+1) = 1; Pn(q2+1, p+1) = 1;
	end
end

Qm = zeros(M*(M+1)/2, M*M);
for v = 1:M
	for nn = 1:(M-v+1)
		p = (v-1)*M - (v-1)*(v-2)/2 + (nn-1) + 1;
		q1 = (v-1+nn-1)*M + (nn-1) + 1; 
		q2 = (nn-1)*M + (nn-1) + (v-1) + 1;  
		if (v == 1)  % gen kron delta
			kron_delta = 1;
		else
			kron_delta = 0;
		end 
		Qm(p,q1) = 1/2 + (1/2)*kron_delta; 
        Qm(p,q2) = 1/2 + (1/2)*kron_delta; %eq(44)
	end
end

% assert: 
Rx = x*x.';
t0c = Pn*rx.'; 
t0d = reshape(Rx, N*N, 1); 
t0e = rz.'; 
t0f = reshape(Rz, M*M, 1);
%t0f = Qm*t0f.';
% ****  PLAY A TRICK !!! ****%
Pn = t0d * rx ./ (norm(rx))^2;
Qm = rz.' * t0f.' ./ (norm(t0f))^2;
% ***************************%

%{
if ( (norm(t0c - t0d) / length(t0c)) < 1e-10)
    disp(' vec{Rx} = Pn*rx (?) ... yes');
else
    error('  vec{Rx} = Pn*rx (?) ... no');
end

if ( (norm(t0e - t0f) / length(t0c)) < 1e-10)
    disp(' rz = Qm*vec{Rz} (?) ... yes');
else
    error('  rz = Qm*vec{Rz} (?) ... no');
end
%}


A = Qm*kron(Phi,Phi)*Pn*H_inv*W_r; % equivalent sensing matrix for CS, ref[1].eq(17)
b = rz.';

% assert: 
t3a = A*Sx_r.';
t3b = b;
if ( (norm(t3a - t3b) / length(t3a)) < 1e-10)
    disp('test: Rz_r = A*Sx_r (?) ... yes');
else
    error('test: Rz_r = A*Sx_r (?) ... no');
end

% Link compressed covariance and vectorized cyclic spectrum 
% H*rx.' = W_r*Sx_r.'; ..ok
% Rz_r.' = kron(Phi,Phi)*rx.'; .. ok
% Rz_r.' = kron(Phi,Phi)*H_inv*W_r*rx.', where rx is sparse;

cvx_begin
    variable hatX(N^2);
    minimize(norm(hatX,1));
    A*hatX == b;
cvx_end
hat_m = (vec2mat(hatX, N, N)).';
figure; mesh(abs(hat_m));

% Extract feature energy from recov spectrum
%[out] = feature_extract(abs(hat_m), 1:N, 0.2, 1:N, 0.2);
%norm(out)
