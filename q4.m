%% Parameters
n = 1000; % number of channel symbols
rho1 = 0.8; % correlation coefficient for state 1
rho2 = 0.2; % correlation coefficient for state 2
pi1 = 0.9; % probability of being in state 1 initially
pi2 = 1 - pi1; % probability of being in state 2 initially
P = [0.95 0.05; 0.05 0.95]; % transition matrix of the Markov chain
SNR_dB = -10:1:10; % range of SNR values in dB

%% Generate transmitted bits and channel symbols
B = randi([0 1], 1, n); % transmitted bits
X = randn(1, n); % uncorrelated channel symbols
S = zeros(1, n); % state sequence of the Markov chain
S(1) = rand(1) < pi1; % initial state
for i = 2:n
    S(i) = rand(1) < P(S(i-1)+1, 2); % generate state transition
end
R = sqrt(rho1)*X + sqrt(1-rho1)*randn(1, n); % channel symbols for state 1
R(S==2) = sqrt(rho2)*X(S==2) + sqrt(1-rho2)*randn(1, sum(S==2)); % channel symbols for state 2

%% Baum-Welch algorithm
[EstTR, EstE] = hmmestimate(R, S, 'Pseudotransitions', 1); % estimate transition and emission probabilities
LogL = hmmdecode(R, EstTR, EstE); % log-likelihood of the observed sequence
for i = 1:numel(SNR_dB)
    SNR = 10^(SNR_dB(i)/10); % convert SNR from dB to linear scale
    var_n = var(R)/SNR; % noise variance
    R_noisy = R + sqrt(var_n)*randn(1, n); % add Gaussian noise to the channel symbols
    LogL_noisy = hmmdecode(R_noisy, EstTR, EstE); % log-likelihood of the noisy sequence
    Pe(i) = sum(abs(B - (LogL_noisy > 0.5)))/n; % bit error rate
end

%% Plot results
figure;
semilogy(SNR_dB, Pe, '-o');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('Bit Error Rate vs SNR');
