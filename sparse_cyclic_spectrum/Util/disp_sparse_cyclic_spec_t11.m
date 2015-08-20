function [Test] = disp_sparse_cyclic_spec_t11(x)
% ***********NOTE**************%
%{
% Header 
clc; clear; close all
addpath('./Util/')
addpath('./Data/')
load gain_attr.mat
disp('load gain_attr: gain.noise gain.sig snr_dB')
load bpsk.mat
x = bpsk(1:64);
%}

load cached_matrix.mat 
%disp('load matrix: Gv_save Dv_save D H W_r H_inv B Pn Qm ')

N = length(x);
%x= gain.sig.*x + gain.noise.*(randn(1,N)); 

x = x.';

% Equivalent cyclic spectrum

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

%generate cyclic spectrum 
Rxc = zeros(N,N);
for v = 0:N-1 
	Gv = Gv_save(:,:,v+1);
	Dv = Dv_save(:,:,v+1);
	Rxc = Rxc + Gv*R*Dv; % ref[1].eq(9)
end
Sx = Rxc*D; %ref[1].eq(10)

% Vectorize xcorr and cyclic spectrum
Sx_r = reshape(Sx, 1, N*N); %reshape the cyclic spectrum

% Compressed sampling the signal
cs.sparse = 16;
cs.ratio = 4;
cs.iter = 100;
cs.N = N;
cs.M = round(cs.N/cs.ratio); % num of sensing points
M = cs.M;
load Phi_16_64.mat
y = Phi*x;
Rz = y*y.';

%vec{Rz}
rz = [];
for v = 0:M-1 
    for n = 0:(M-1-v)
        rz = [rz, y(1+n)*y(1+n+v)];
	end
end

A = Qm*kron(Phi,Phi)*Pn*H_inv*W_r; % equivalent sensing matrix for CS, ref[1].eq(17)
b = rz.';

% assert: 
t3a = A*Sx_r.';
t3b = b;

% Link compressed covariance and vectorized cyclic spectrum 
% H*rx.' = W_r*Sx_r.'; ..ok
% rz.' = Qm*kron(Phi,Phi)*Pn*rx.'; .. ok

lambda_opt = 0.05;
cvx_begin quiet
	variable hatX(N^2);
	minimize(lambda_opt.*norm(hatX,1) + norm(b-A*hatX)); %re. eq.21
cvx_end
hat_m = (vec2mat(hatX, N, N)).';

%CS based noise covariance
load feature_matrix.mat %feature_mask Sigma J

%Test-Statistic
c = J*hatX;   % abs(vec2mat(c, N, N)).' = feature_mask.*abs(hat_m); .. ok
Test = c'*Sigma*c;  
