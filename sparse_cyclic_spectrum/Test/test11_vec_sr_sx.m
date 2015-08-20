% ***********NOTE**************%
%function [hat_m_M, Sx_M, feature_mask] = disp_sparse_cyclic_spec(x)
clc; clear; close all
addpath('./Util/')
addpath('./Data/')

% Header 
load gain_attr.mat
disp('load gain_attr: gain.noise gain.sig snr_dB')
load bpsk.mat
x = bpsk(1:64);
load cached_matrix.mat 
disp('load matrix: Gv_save Dv_save D H W_r H_inv B Pn Qm ')

N = length(x);
x= gain.sig.*x + gain.noise.*(randn(1,N)); 

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
%{
Phi = pn_gen(cs.M,cs.N);
Phi_mask = zeros(cs.M, cs.N);
for i = 1:M
    if i*cs.ratio > N
        Phi_mask(i, ((i-1)*cs.ratio+1) : N) = ones(1, N- (i-1)*cs.ratio );
    else    
        Phi_mask(i, ((i-1)*cs.ratio+1) : (i*cs.ratio)) = ones(1, cs.ratio);
    end    
end
Phi = Phi.*Phi_mask;
%}
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

if (norm(imag(t3a)-imag(t3b)) < 1e-10) && (norm(real(t3a)-real(t3b)) < 1e-10)
    disp('test: Rz_r = A*Sx_r (?) ... yes');
else
    error('test: Rz_r = A*Sx_r (?) ... no');
end

% Link compressed covariance and vectorized cyclic spectrum 
% H*rx.' = W_r*Sx_r.'; ..ok
% rz.' = Qm*kron(Phi,Phi)*Pn*rx.'; .. ok

%[hatX, ~] = cosamp(b, A, cs.sparse, cs.iter);

% test lambda_opt
%{
lambda_i = 1;
for lambda = 0.01: 0.02: 0.1
	cvx_begin
		variable hatX(N^2);
		minimize(lambda.*norm(hatX,1) + norm(b-A*hatX));
	cvx_end
	hatX_R(:,lambda_i) = hatX;
	lambda_i = lambda_i + 1;
    l1norm(lambda_i) = norm(hatX,1);
    l2norm(lambda_i) = norm(A*hatX-b);
end
figure; plot( l1norm, l2norm, 'o' );
xlabel( 'norm(x,1)' );
ylabel( 'norm(A*x-b)' );
grid on
%}
lambda_opt = 0.05;
cvx_begin
	variable hatX(N^2);
	minimize(lambda_opt.*norm(hatX,1) + norm(b-A*hatX)); %re. eq.21
cvx_end
hat_m = (vec2mat(hatX, N, N)).';

% Extract feature energy from recov spectrum
feature_mask = zeros(N,N);
for i = 0:N-1
    for j = 0:N-1
        if (j <= 2*i+14) && (j >= 2*(i-11))  
            feature_mask(j+1,i+1) = 1;
        end
        if (j >= -2*(i-53)) && (j <= -2*(i-63) + 15)
            feature_mask(j+1,i+1) = 1;
        end
    end
end

% Show
figure; mesh(abs(Sx));
figure; mesh(abs(hat_m));
figure; mesh(feature_mask.*abs(Sx));
figure; mesh(feature_mask.*abs(hat_m));

%CS based noise covariance
L = 10; %test group number
for l = 1:L
	x_test = gain.sig.*x + gain.noise.*(randn(N,1)); 
	y_test = Phi*x_test;
	Rz_test = y_test * y_test';
	rz_test(:,l) = Qm*reshape(Rz_test,M*M,1);
end
rz_test_mean = (1/L).*sum(rz_test,2); 
Sigma_z = zeros(length(rz_test_mean), length(rz_test_mean));
for l = 1:L
	Sigma_z = Sigma_z + (rz_test(:,l) - rz_test_mean) * (rz_test(:,l) - rz_test_mean)'; %eq.45
end
Sigma_z = (1/L).*Sigma_z; 
T = ((A'*A+lambda_opt.*eye(N^2))^(-1)) * A';
J_tmp = reshape(feature_mask, 1, N*N);
J = diag(J_tmp);
%Sigma
Sigma = J*T*Sigma_z*T'*J'; %eq.47

%Test-Statistic
c = J*hatX;   % abs(vec2mat(c, N, N)).' = feature_mask.*abs(hat_m); .. ok
Test = c'*Sigma*c 
save feature_matrix.mat feature_mask Sigma J