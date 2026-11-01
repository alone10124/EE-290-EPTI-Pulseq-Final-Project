function [dB, fit_error, phase] = dB_fitting_JumpCorrect(phase, TEs, mask, unwrap_flag)
% Fitting of the B0 inhomogeneity field (dB) using multi-echo phase data
% Fuyixue Wang, 2018

% updated, clean up
% Zijing Dong, 2019

% update, vectorized
% Daniel Haenelt, 2026

if unwrap_flag == 1
    phase = unwrap(phase, [], 3);
end

% Normalize phase
phase = phase ./ (2*pi);

[nx, ny, nt] = size(phase);

% Design matrix
A = [ones(nt, 1), TEs(1:nt)];

% Reshape data
phase2 = reshape(phase, [], nt).';  % [nt x nx*ny]
mask1d = mask(:).';

% Select masked voxels
Y = double(phase2(:, mask1d));  % [nt x masked voxels]

% Least-squares fit
param = A \ Y;

% Fit signal and error
Yfit = A * param;
err = Yfit - Y;
rmse = sqrt(mean(err.^2, 1));

% Allocate outputs
dB_map = zeros(1, nx*ny);
fit_error_map = zeros(1, nx*ny);

dB_map(mask1d) = param(2, :);
fit_error_map(mask1d) = rmse;

% Reshape back to image space
dB = reshape(dB_map, nx, ny);
dB(isnan(dB)) = 0;
fit_error = reshape(fit_error_map, nx, ny);

% Restore phase to radians
phase = phase * (2*pi);

end
