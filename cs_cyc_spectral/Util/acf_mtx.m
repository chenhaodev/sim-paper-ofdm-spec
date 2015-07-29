function [rxx] = acf_mtx(x)
%autocorrelation matrix generate for a vector N * 1
% https://en.wikipedia.org/wiki/Autocorrelation_matrix
% http://www.mathworks.com/matlabcentral/answers/49172-autocorrelation-matrix-from-a-vector

N = length(x);

[xc,~] = xcorr(x,x,length(x)-1,'biased');
r = xc(length(x):end);
rxx = toeplitz(r,conj(r)); 
%{
d_alpha = fs/N; % freq resolution
alpha = 0:d_alpha:fs-d_alpha; % cyclic resolution
a_len = length(alpha); 

f_len = N; 
f = -(fs/2-d_alpha*floor(N/2)) + d_alpha*(0:N-1); % freq sample location

Rxx = fft(rxx);
for i = 1:N
    RXX(:,i) = simple_fft(Rxx(:,i), a_len, d_alpha);
end

Spec = abs(RXX);
mesh(f, alpha, Spec); 
axis tight;
xlabel('f'); ylabel('a'); 
%}