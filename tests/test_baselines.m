%% Test Baseline Algorithms: ZF and SQP
% Quick test to verify the baseline algorithms work correctly

clear; close all; clc;

%% Setup
fprintf('=== Testing Baseline Algorithms ===\n\n');

% Add paths
addpath(genpath('src'));

% System parameters
num_antennas = 8;
num_users = 4;
num_realizations = 5;

% Set seed for reproducibility
set_global_seed(42);

%% Generate parameters
config.noise_power = ones(num_users, 1);
config.gamma = ones(num_users, 1) * 2.0;  % 2.0 linear = 3 dB

% SQP-specific parameters (optional, using defaults if not specified)
config.sqp_max_iter = 500;
config.sqp_tol = 1e-6;
config.sqp_display = 'off';  % Change to 'iter' to see optimization progress

%% Run tests
fprintf('Configuration:\n');
fprintf('  Antennas: %d, Users: %d\n', num_antennas, num_users);
fprintf('  QoS target: %.1f (%.1f dB)\n', config.gamma(1), 10*log10(config.gamma(1)));
fprintf('  Noise power: %.1f\n\n', config.noise_power(1));

zf_powers = zeros(num_realizations, 1);
sqp_powers = zeros(num_realizations, 1);
zf_times = zeros(num_realizations, 1);
sqp_times = zeros(num_realizations, 1);

for r = 1:num_realizations
    fprintf('--- Realization %d/%d ---\n', r, num_realizations);
    
    % Generate random Rayleigh fading channel
    H = (randn(num_antennas, num_users) + 1j*randn(num_antennas, num_users)) / sqrt(2);
    
    %% Test Zero-Forcing
    [W_zf, metrics_zf] = zf_beamformer(H, config);
    
    fprintf('ZF Beamformer:\n');
    fprintf('  Power: %.4f W (%.2f dB)\n', metrics_zf.final_power, metrics_zf.power_db);
    fprintf('  Min SNR: %.4f (%.2f dB)\n', metrics_zf.min_snr, metrics_zf.min_snr_db);
    fprintf('  Feasible: %d\n', metrics_zf.feasible);
    fprintf('  Time: %.6f s\n', metrics_zf.solve_time);
    fprintf('  Status: %s\n', metrics_zf.status_message);
    
    zf_powers(r) = metrics_zf.final_power;
    zf_times(r) = metrics_zf.solve_time;
    
    %% Test SQP
    [W_sqp, metrics_sqp] = sqp_beamformer(H, config);
    
    fprintf('SQP Beamformer:\n');
    fprintf('  Power: %.4f W (%.2f dB)\n', metrics_sqp.final_power, metrics_sqp.power_db);
    fprintf('  Min SNR: %.4f (%.2f dB)\n', metrics_sqp.min_snr, metrics_sqp.min_snr_db);
    fprintf('  Feasible: %d\n', metrics_sqp.feasible);
    fprintf('  Converged: %d\n', metrics_sqp.converged);
    fprintf('  Iterations: %d\n', metrics_sqp.iterations);
    fprintf('  Time: %.6f s\n', metrics_sqp.solve_time);
    fprintf('  Status: %s\n', metrics_sqp.status_message);
    
    sqp_powers(r) = metrics_sqp.final_power;
    sqp_times(r) = metrics_sqp.solve_time;
    
    fprintf('\n');
end

%% Summary statistics
fprintf('=== Summary Statistics (%d realizations) ===\n\n', num_realizations);

fprintf('Zero-Forcing:\n');
fprintf('  Average power: %.4f W\n', mean(zf_powers));
fprintf('  Average time: %.6f s\n', mean(zf_times));

fprintf('\nSQP:\n');
fprintf('  Average power: %.4f W\n', mean(sqp_powers));
fprintf('  Average time: %.6f s\n', mean(sqp_times));

fprintf('\nSQP vs ZF:\n');
fprintf('  Power ratio: %.2fx (SQP should be <= ZF)\n', mean(sqp_powers) / mean(zf_powers));
fprintf('  Time ratio: %.2fx (SQP is typically slower)\n', mean(sqp_times) / mean(zf_times));

%% Validation
fprintf('\n=== Validation ===\n');
if mean(sqp_powers) <= mean(zf_powers) * 1.05  % Allow 5% tolerance
    fprintf('✓ SQP power <= ZF power (as expected for optimal vs suboptimal)\n');
else
    fprintf('✗ Warning: SQP power > ZF power (unexpected)\n');
end

if all(isfinite(zf_powers)) && all(isfinite(sqp_powers))
    fprintf('✓ All solutions are finite\n');
else
    fprintf('✗ Warning: Some solutions are infinite/NaN\n');
end

fprintf('\nTest completed successfully!\n');
