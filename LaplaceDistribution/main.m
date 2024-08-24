% Number of samples
N = 10000;

% Correlation matrix for two Gaussian variables
rho = 0.5; % Correlation coefficient
Sigma = [1, rho; rho, 1];

% Generate correlated Gaussian random variables
R = chol(Sigma);
gaussSamples = repmat([0,0], N, 1) + randn(N,2) * R;

% Transform Gaussian to Laplacian
% Using the quantile transformation method
laplaceSamples(:,1) = sign(gaussSamples(:,1)) .* log(1 - 2 * abs(normcdf(gaussSamples(:,1)) - 0.5));
laplaceSamples(:,2) = sign(gaussSamples(:,2)) .* log(1 - 2 * abs(normcdf(gaussSamples(:,2)) - 0.5));

% Visualization with a scatter plot
scatter(laplaceSamples(:,1), laplaceSamples(:,2), 10, 'filled');
xlabel('Variable 1');
ylabel('Variable 2');
title('Scatter Plot of Correlated Laplacian Random Variables');
axis equal;
grid on;

% Assuming laplaceSamples is obtained from the previous step

% Estimate the density
n = 100;
[xi, yi] = meshgrid(linspace(min(laplaceSamples(:,1)), max(laplaceSamples(:,1)), n), ...
                    linspace(min(laplaceSamples(:,2)), max(laplaceSamples(:,2)), n));
f = ksdensity(laplaceSamples, [xi(:) yi(:)], 'Kernel', 'normal', 'Bandwidth', 0.5);
f = reshape(f, n, n);

% Contour plot of the estimated density
contour(xi, yi, f, 20); % Adjust the number of contours as needed
xlabel('Variable 1');
ylabel('Variable 2');
title('Contour Plot of Correlated Laplacian Random Variables');
colorbar; % Optional: Adds a colorbar to indicate density levels
