clc
clear all;
N=10000;%独立随机数实现的数目
OR=20;%观察频率（Hz）
M=256;%多普勒滤波器的阶数
Dop_res=0.1;%SUI参数中的多普勒判决（Hz）（在重复采样进程中）
res_accu=20;%重复采样进程的精确度
%%%%%SUI信道参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=[0 -5 -10];%每一阶的功率衰减
K=[1 0 0];%K因子
tau=[0.0 0.4 0.9];%每个抽头的时延
Dop=[0.4 0.3 0.5];%最大多普勒频移参数
ant_corr=0.4;%天线相关性
Fnorm=-1.5113;%增益归一化因子
%%%%%%%%%%%计算m的值%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 10.^(P/10);%计算线性功率
s2 = P./(K+1); % 计算方差
m2 = P.*(K./(K+1)); % 计算常数功率
m = sqrt(m2); % 计算常数部分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmsdel = sqrt( sum(P.*(tau.^2))/sum(P) - (sum(P.*tau)/sum(P))^2 );
fprintf('rms delay spread %6.3f μs\n', rmsdel);
%%%%%%%%%在特定功率下计算伦琴信道系数%%%%%%%%%%%%%%%%
L = length(P); % 阶数
paths_r = sqrt(1/2)*(randn(L,N) + j*randn(L,N)).*((sqrt(s2))' * ones(1,N)); %L*N矩阵每阶的数据噪声
paths_c = m' * ones(1,N);%常数部分

for p = 1:L 
    D = Dop(p) / max(Dop) / 2; % 归一化最大多普勒频移 
    f0 = [0:M*D]/(M*D); % 频率因子 
    PSD = 0.785*f0.^4 - 1.72*f0.^2 + 1.0; % PSD估计 
    filt = [ PSD(1:end-1) zeros(1,M-2*M*D) PSD(end:-1:2) ]; % S(f) 
    filt = sqrt(filt); %从S(f)到|H(f)| 
    filt = ifftshift(ifft(filt)); % 获得脉冲响应 
    filt = real(filt); % 寻找实数滤波器
    filt = filt / sqrt(sum(filt.^2)); %归一化滤波器
    path = fftfilt(filt, [ paths_r(p,:) zeros(1,M) ]); 
    paths_r(p,:) = path(1+M/2:end-M/2); 
end; 
paths = paths_r + paths_c;%路径输出数据
Pest = mean(abs(paths).^2, 2);
fprintf('tap mean power level: %0.2f dB\n', 10*log10(Pest));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SR = max(Dop)*2; % 精确采样率
m = lcm(SR/Dop_res, OR/Dop_res); %求最小公倍数
P = m/SR*Dop_res; % 分子
Q = m/OR*Dop_res; % 分母 
paths_OR = zeros(L,ceil(N*P/Q)); % 创造新矩阵
for p=1:L
    paths_OR(p,:) = resample(paths(p,:), P, Q, res_accu);
end;
Pest = mean(abs(paths_OR).^2, 2);
fprintf('tap mean power level: %0.2f dB\n', 10*log10(Pest));
paths_OR1=10*log10(paths_OR);
NN=length(paths_OR(1,:));
y1=abs(paths_OR).^2;
y2=10*log10(y1);%转换为dB
t=60;%时间长度
x=1:NN;
y=x/OR;
plot(y,y2(1,:),y,y2(2,:),'--',y,y2(3,:),'-.')
axis([1,60,-60,10])
legend('tab1','tab2','tab3')
xlabel('time(s), Fs = 20Hz')
ylabel('Amplitude dB')


