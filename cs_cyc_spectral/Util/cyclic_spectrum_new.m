function [Spec, f, alpha] = cyclic_spectrum_new(x, N, fs, M, opt)
% cyclic spectrum analysis
% x: signal 
% N: samples <= len(x) 
% fs: sample rate, [-fs/2 ~ fs/2]
% M: window length 
% author: chenhaomails@gmail.com

win = 'hamming';

d_alpha = fs/N; % freq resolution
alpha = 0:d_alpha:fs-d_alpha; % cyclic resolution
a_len = length(alpha); 

%f_len = floor(N/M-1)+1; 
f_len = ceil(N/M-1)+1; 
f = -(fs/2-d_alpha*floor(M/2)) + d_alpha*M*(0:f_len-1); % freq sample location

S = zeros(a_len, f_len); 
i = 1; 

% signal fft
X = fftshift(fft(x(1:N))); 
X = X';

%% Loop
for alfa = alpha

    interval_f_N = round(alfa/d_alpha);
    f_N = floor((N-interval_f_N-M)/M)+1; % window num ~= N/M
    f_N_m = ceil((N-interval_f_N-M)/M)+1;
    add_zero_valid = 0;
    add_zero_num = 0;
    if f_N < f_N_m
       add_zero_valid = 1;
       add_zero_num = M* (f_N_m - f_N);
    end
    
    % window generate
    g = feval(win, M); % return an M-point window  
    window_M = g(:, ones(f_N,1));
    t = 1:M*f_N;
    t = reshape(t, M, f_N);

    % spectral correlation
    X1 = X(t);%.*window_M; 
    X2 = X(t+interval_f_N);%.*window_M; 
    if add_zero_valid ==1
        X1 = [X1 zeros(add_zero_num,1)];
        X2 = [X2 zeros(add_zero_num,1)];
    end    
    %St = conj(X1).*X2;  
    St = conj(X2).*X1;  
    St = mean(St, 1); % T average
    %S(i, floor((f_len-f_N)/2)+(1:f_N)) = St/N; %move St to central
    S(i, floor((f_len-f_N_m)/2)+(1:f_N_m)) = St/N; %move St to central
    i = i+1;

    
end

Spec = abs(S);

% figure
if strcmpi(opt,'show')
    mesh(f, alpha, Spec); 
    axis tight;
    xlabel('f'); ylabel('a');    
end
