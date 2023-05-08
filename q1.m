% Set parameters
n = 10000; % number of channel symbols
rho_values = -1:0.1:1; % values of rho to be tested
M = 100; % number of Monte Carlo trials for each rho value

% Initialize error rate matrix
Pe = zeros(length(rho_values), 1);

% Loop over rho values
for i = 1:length(rho_values)
    rho = rho_values(i);
    for j = 1:M
        % Generate random bits
        b = randi([0, 1], 1, n);

        % Generate random noise with correlation rho
        w = sqrt(1 - rho^2) * randn(1, n) + rho * randn(1, n);

        % Generate received signal
        y = b + w;

        % Make detection based on correlation coefficient
        if rho == 0
            % Check first component of y
            hat_b = (y(1) > 0);
        else
            % Check sign of y
            hat_b = (sum(sign(y)) > 0);
        end

        % Calculate error rate
        Pe(i) = Pe(i) + (hat_b ~= b(1));
    end

    % Average error rate over Monte Carlo trials
    Pe(i) = Pe(i) / M;
end

% Plot error rate as a function of rho
figure;
plot(rho_values, Pe, 'o-');
xlabel('Correlation coefficient \rho');
ylabel('Error rate P_e');
ylim([0, 0.5]);
