clear;
close all;
clc;

% Number of samples per realization
N = 1024;

% Number of realizations
numRealizations = 1; % for example, generate 10 realizations

% Time vector (assuming a sampling interval of 1)
t = (0:N-1)';

% Discretize the autocorrelation function R(tau) = exp(-abs(tau))
l = 100;
tau = t - N/2;
R = exp(-((tau.^2)/(2 * l^2)));
R = tau;

% Compute the Power Spectral Density (PSD) using FFT
PSD = fftshift(fft(ifftshift(R)));

% Initialize matrix to store all realizations
realizations = zeros(N, numRealizations);

% Generate multiple realizations
for k = 1:numRealizations
    % Generate random phase
    randomPhase = exp(2i * pi * rand(size(PSD)));

    % Generate the signal in the frequency domain by combining PSD and random phase
    S_freq = sqrt(PSD) .* randomPhase;

    % Convert back to the time domain to get one realization
    signal = real(ifft(ifftshift(S_freq)));
    
    % Store the realization
    realizations(:, k) = signal;
end

% Plot the generated stochastic process realizations
figure;
plot(t, realizations);
xlabel('t(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex');
xlim([0,1000])
%append the value of l to the title
title(strcat('Ensemble of Stochastic Process Realizations, $l = ', num2str(l), "$") , 'Interpreter','latex');
exportgraphics(gcf, 'stochastic_process_realizations.png')
