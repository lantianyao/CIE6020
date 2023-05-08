clear all; close all; clc;

% ��������
N = 10000; % ���������
SNRdB = 0:2:16; % ����ȷ�Χ
sigmas = 10.^(-SNRdB/20); % ��׼��
rho = 0.5; % δ֪����rho
m = 50; % ѵ�����г���

% ����ѵ������
trainSeq = randi([0 1], 1, m);

% ��ʼ�������ʾ���
Pe = zeros(length(sigmas), 1);

for i = 1:length(sigmas)
    sigma = sigmas(i);
    fprintf('Running simulation %d of %d for sigma = %.4f\n', i, length(sigmas), sigma);
    
    % ���ɷ��ͱ���
    bits = randi([0 1], 1, N);

    % QPSK����
    modSig = 1 - 2*bits(1:2:end) + 1i*(1 - 2*bits(2:2:end));
    
    % ��Ӹ�˹������
    rxSig = modSig + sigma*(randn(size(modSig)) + 1i*randn(size(modSig)));
    
    % ����rho
    r = sum(rxSig(1:m) .* conj(rxSig(2:m+1)));
    estRho = abs(r)/m;

    % �����
    bitsHat = zeros(1, N);
    for j = m+1:N/2
        xHat = rxSig(j)/sqrt(estRho);
        bitsHat(2*j-1) = real(xHat) < 0;
        bitsHat(2*j) = imag(xHat) < 0;
    end
    
    % ����������
    Pe(i) = sum(bits ~= bitsHat)/N;
end

% ����Pe��SNR�Ĺ�ϵͼ
figure;
semilogy(SNRdB, Pe, '-o');
grid on;
xlabel('SNR (dB)');
ylabel('Pe');
title(['Pe vs. SNR for QPSK with training sequence (m = ', num2str(m), ')']);

