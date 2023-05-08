clear all; close all; clc;

% 参数设置
N = 10000; % 传输比特数
SNRdB = 0:2:16; % 信噪比范围
sigmas = 10.^(-SNRdB/20); % 标准差
rho = 0.5; % 未知参数rho
m = 50; % 训练序列长度

% 生成训练序列
trainSeq = randi([0 1], 1, m);

% 初始化误码率矩阵
Pe = zeros(length(sigmas), 1);

for i = 1:length(sigmas)
    sigma = sigmas(i);
    fprintf('Running simulation %d of %d for sigma = %.4f\n', i, length(sigmas), sigma);
    
    % 生成发送比特
    bits = randi([0 1], 1, N);

    % QPSK调制
    modSig = 1 - 2*bits(1:2:end) + 1i*(1 - 2*bits(2:2:end));
    
    % 添加高斯白噪声
    rxSig = modSig + sigma*(randn(size(modSig)) + 1i*randn(size(modSig)));
    
    % 估计rho
    r = sum(rxSig(1:m) .* conj(rxSig(2:m+1)));
    estRho = abs(r)/m;

    % 解调器
    bitsHat = zeros(1, N);
    for j = m+1:N/2
        xHat = rxSig(j)/sqrt(estRho);
        bitsHat(2*j-1) = real(xHat) < 0;
        bitsHat(2*j) = imag(xHat) < 0;
    end
    
    % 计算误码率
    Pe(i) = sum(bits ~= bitsHat)/N;
end

% 绘制Pe与SNR的关系图
figure;
semilogy(SNRdB, Pe, '-o');
grid on;
xlabel('SNR (dB)');
ylabel('Pe');
title(['Pe vs. SNR for QPSK with training sequence (m = ', num2str(m), ')']);

