clc; clear; close all

addpath('./Util/')
addpath('./Data/')


% Header 

sig.type = 'bpsk'; % 'fsk'
sig.fs = 1;
sig.M = 1;

if strcmpi(sig.type,'fsk') % default signal
	load fsk.mat
    sig.x = fsk_real(1:64);
elseif strcmpi(sig.type,'bpsk') % default signal
	load bpsk.mat
    sig.x = bpsk_w(1:64);
else
	error('signal type not exist!!');
end

sig.x= randn(1,64); %!!!!!!

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

%rx generate, ref[1].eq(8)
rx = [];
for v = 0:N-1 
    for n = 0:(N-1-v)
        rx = [rx, x(1+n)*x(1+n+v)];
	end
end

%R generation, ref[1].eq(7)
R = zeros(N,N);
ptr = 1;
for i = 0:N-1
    num = N - i;
    R(:,i+1) = [rx(ptr:ptr+num-1) zeros(1,i)]';
    ptr = ptr + num;
end

%vec{R} generate, ref[1].eq(8)
B = [];
tmp = eye(N*N);
for v = 0:N-1 
    for n = 0:(N-1-v)
		q = v*N + n;
		B = [B, tmp(:, q+1)];
	end
end

%{
% assert: 
t0c = B*rx.'; 
t0d = reshape(R, N*N, 1);
if ( (norm(t0c - t0d) / length(t0c)) < 1e-10)
    disp(' vec{R} = B*rx (?) ... yes');
else
    error(' vec{R} = B*rx (?) ... no');
end
%}

%generate cyclic spectrum 
Rxc = zeros(N,N);
for v = 1:N 
	Gv = zeros(N,N);
	Dv = zeros(N,N);
	for nn = 1:N
		for alph = 1:N
			%Gv(nn, alph) = exp(-2*pi*sqrt(-1)*(alph-1)*(nn-1+((v-1)/2))/N); % recall dftmtx(N): each W(j,k)= exp(-2*pi*sqrt(-1)*(j-1)*(k-1)/N); 
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
t1a = W_r*Sx_r.'; 
t1b = H*rx.';
if ( (norm(t1a - t1b) / length(t1a)) < 1e-10)
    disp('test: H*rx == W_r*Sx_r (?) ... yes');
else
    error('test: H*rx == W_r*Sx_r (?) ... no');
end

% Compressed sampling the signal
cs.sparse = 16;
cs.ratio = 4;
cs.iter = 32;
cs.N = N;
cs.M = round(cs.N/cs.ratio); % num of sensing points
M = cs.M;
Phi = randn(cs.M, cs.N);
y = Phi*x;
Rz = y*y.';

%vec{Rz}
rz = [];
for v = 0:M-1 
    for n = 0:(M-1-v)
        rz = [rz, y(1+n)*y(1+n+v)];
	end
end

%gen Pn, Qm
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
t2a = Pn*rx.'; 
t2b = reshape(Rx, N*N, 1); 
if ( (norm(t2a - t2b) / length(t2a)) < 1e-10)
    disp(' vec{Rx} = Pn*rx (?) ... yes');
else
    error('  vec{Rx} = Pn*rx (?) ... no');
end
t2c = rz.'; 
t2d = Qm*reshape(Rz, M*M, 1);
if ( (norm(t2c - t2d) / length(t2c)) < 1e-10)
    disp(' rz = Qm*vec{Rz} (?) ... yes');
else
    error('  rz = Qm*vec{Rz} (?) ... no');
end

% Link compressed covariance and vectorized cyclic spectrum 
% H*rx.' = W_r*Sx_r.'; ..ok
% rz.' = Qm*kron(Phi,Phi)*Pn*rx.'; .. ok
% Hence, equivalent sensing matrix for CS, ref[1].eq(17):

% assert: 
t3a = Qm*kron(Phi,Phi)*Pn*rx.';
t3b = rz.';
if ( (norm(t3a - t3b) / length(t3a)) < 1e-10)
    disp('test: rz = Qm*kron(Phi,Phi)*Pn*rx (?) ... yes');
else
    error('test: rz = Qm*kron(Phi,Phi)*Pn*rx (?) ... no');
end
A = Qm*kron(Phi,Phi)*Pn*H_inv*W_r;
b = rz.';

    cvx_begin
        variable x(N^2);
        minimize( norm(A*x-b)+5.4556*norm(x,1) );
    cvx_end
hat_m = (vec2mat(x, N, N)).';
figure; mesh(abs(hat_m));


%{
cvx_begin
    variable hatX(N^2);
    minimize(norm(hatX,1));
    A*hatX == b;
cvx_end
hat_m = (vec2mat(hatX, N, N)).';
figure; mesh(abs(hat_m));
%figure; mesh(abs(vec2mat(W_r*hatX, N,N)*D));
%}

% Extract feature energy from recov spectrum
%[out] = feature_extract(abs(hat_m), 1:N, 0.2, 1:N, 0.2);
%norm(out)

%{
gamma = logspace( -2, 2, 20 );
l2norm = zeros(size(gamma));
l1norm = zeros(size(gamma));
fprintf( 1, '   gamma       norm(x,1)    norm(A*x-b)\n' );
fprintf( 1, '---------------------------------------\n' );
for k = 1:length(gamma),
    fprintf( 1, '%8.4e', gamma(k) );
    cvx_begin
        variable x(N^2);
        minimize( norm(A*x-b)+gamma(k)*norm(x,1) );
    cvx_end
    l1norm(k) = norm(x,1);
    l2norm(k) = norm(A*x-b);
    fprintf( 1, '   %8.4e   %8.4e\n', l1norm(k), l2norm(k) );
end
plot( l1norm, l2norm );
xlabel( 'norm(x,1)' );
ylabel( 'norm(A*x-b)' );
grid on
%}