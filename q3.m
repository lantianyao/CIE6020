%% Part (c): EM algorithm to estimate rho
clear all; close all; clc;

%% Parameters
n = 10000; % number of symbols
M = 100; % number of Monte Carlo trials
SNRdBs = -5:2:15; % range of SNR values
SNRs = 10.^(SNRdBs/10);
rho_true = 0.5; % true value of rho

%% Initialization
Pe_EM = zeros(length(SNRs), 1);
Pe_TS = zeros(length(SNRs), 1);

%% Main loop
for i = 1:length(SNRs)
    SNR = SNRs(i);
    sigma2 = 0.5 ./ SNR;
    
    for j = 1:M
        %% Generate channel and symbols
        h = sqrt(rho_true)*randn(1, n) + sqrt(1 - rho_true)*randn(1, n);
        s = sign(randn(1, n));
        x = sqrt(sigma2).*randn(1, n);
        r = h.*s + x;
        
        %% EM algorithm
        rho_EM = 0.5; % initial guess
        for k = 1:100 % EM iterations
            r_hat = sqrt(rho_EM).*s + sqrt(1 - rho_EM).*x;
            z = r.*conj(r_hat) ./ (abs(r_hat).^2);
            rho_EM = mean(z);
        end
        
        %% Detection with training sequence
        m = 100; % length of training sequence
        r_TS = r(m+1:end);
        s_TS = s(m+1:end);
        z_TS = r_TS.*conj(s_TS) ./ (abs(s_TS).^2);
        rho_TS = mean(z_TS);
        
        %% Calculate Pe
        if rho_EM >= 0
            Pe_EM(i) = Pe_EM(i) + (1 - qfunc(sqrt(SNR/(1 + SNR)*abs(rho_EM))));
        else
            Pe_EM(i) = Pe_EM(i) + (1 - qfunc(sqrt(SNR/(1 + SNR)*abs(rho_EM))) - exp(-SNR/(1 + SNR)*(1 - abs(rho_EM))));
        end
        
        if rho_TS >= 0
            Pe_TS(i) = Pe_TS(i) + (1 - qfunc(sqrt(SNR/(1 + SNR)*abs(rho_TS))));
        else
            Pe_TS(i) = Pe_TS(i) + (1 - qfunc(sqrt(SNR/(1 + SNR)*abs(rho_TS))) - exp(-SNR/(1 + SNR)*(1 - abs(rho_TS))));
        end
    end
    
    %% Average Pe over Monte Carlo trials
    Pe_EM(i) = Pe_EM(i) / M;
    Pe_TS(i) = Pe_TS(i) / M;
end

%% Plot Pe versus SNR
figure;
semilogy(SNRdBs, Pe_EM, '-o', 'LineWidth', 2);
hold on;
semilogy(SNRdBs, Pe_TS, '-^', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Pe');
legend('EM algorithm', 'Training sequence');
