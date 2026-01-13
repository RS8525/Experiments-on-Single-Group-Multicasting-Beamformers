% plot_example_results.m
% Example plotting script for multicast beamforming results
%
% This script demonstrates:
% - Loading saved results
% - Generating performance figures
% - Saving figures to results/figures/

% Add paths
addpath(genpath('../src'));

% Clear workspace
clear; close all; clc;

%% Load results
fprintf('Loading results...\n');

% Specify result file (modify this path to your actual result file)
result_file = '../results/raw/example_experiment_20231215_143022.mat';

% If the file doesn't exist, try to find the most recent one
if ~isfile(result_file)
    result_dir = '../results/raw';
    files = dir(fullfile(result_dir, 'example_experiment_*.mat'));
    
    if isempty(files)
        error('No result files found in %s', result_dir);
    end
    
    % Sort by date and get most recent
    [~, idx] = sort([files.datenum], 'descend');
    result_file = fullfile(result_dir, files(idx(1)).name);
    fprintf('Using most recent result file: %s\n', files(idx(1)).name);
end

load(result_file, 'results');

%% Extract data
snr_db = 10 * log10(results.snr_ff_c2);
rate = results.rate_ff_c2;
iterations = results.iterations_ff_c2;

K = size(snr_db, 1);  % Number of users
num_realizations = size(snr_db, 2);

%% Create output directory
output_dir = '../results/figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Figure 1: Sum Rate Distribution
figure('Position', [100, 100, 800, 600]);

subplot(2, 2, 1);
histogram(rate, 20, 'Normalization', 'probability');
xlabel('Sum Rate (bits/s/Hz)');
ylabel('Probability');
title('Sum Rate Distribution');
grid on;

subplot(2, 2, 2);
plot(rate, 'b.-', 'LineWidth', 1.5);
xlabel('Realization');
ylabel('Sum Rate (bits/s/Hz)');
title('Sum Rate vs Realization');
grid on;

subplot(2, 2, 3);
boxplot(snr_db');
xlabel('User Index');
ylabel('SNR (dB)');
title('SNR Distribution per User');
grid on;

subplot(2, 2, 4);
histogram(iterations, 20);
xlabel('Iterations');
ylabel('Count');
title('Convergence Iterations');
grid on;

sgtitle(sprintf('FF-C2 Performance (M=%d, K=%d)', ...
        results.config.system.num_antennas, ...
        results.config.system.num_users));

% Save figure
fig_file = fullfile(output_dir, 'example_performance_summary.png');
saveas(gcf, fig_file);
fprintf('Figure saved: %s\n', fig_file);

%% Figure 2: SNR per User
figure('Position', [150, 150, 800, 500]);

% Mean SNR per user
mean_snr = mean(snr_db, 2);
std_snr = std(snr_db, 0, 2);

errorbar(1:K, mean_snr, std_snr, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('User Index');
ylabel('SNR (dB)');
title('Mean SNR per User with Standard Deviation');
grid on;
ylim([min(mean_snr - std_snr) - 2, max(mean_snr + std_snr) + 2]);

% Save figure
fig_file = fullfile(output_dir, 'example_snr_per_user.png');
saveas(gcf, fig_file);
fprintf('Figure saved: %s\n', fig_file);

%% Figure 3: CDF of Sum Rate
figure('Position', [200, 200, 600, 500]);

[f, x] = ecdf(rate);
plot(x, f, 'b-', 'LineWidth', 2);
xlabel('Sum Rate (bits/s/Hz)');
ylabel('CDF');
title('Cumulative Distribution of Sum Rate');
grid on;

% Add median line
median_rate = median(rate);
hold on;
plot([median_rate, median_rate], [0, 1], 'r--', 'LineWidth', 1.5);
legend('CDF', sprintf('Median (%.3f)', median_rate), 'Location', 'southeast');

% Save figure
fig_file = fullfile(output_dir, 'example_rate_cdf.png');
saveas(gcf, fig_file);
fprintf('Figure saved: %s\n', fig_file);

%% Display statistics
fprintf('\n========== Statistics Summary ==========\n');
fprintf('Sum Rate:\n');
fprintf('  Mean:   %.4f bits/s/Hz\n', mean(rate));
fprintf('  Median: %.4f bits/s/Hz\n', median(rate));
fprintf('  Std:    %.4f bits/s/Hz\n', std(rate));
fprintf('  Min:    %.4f bits/s/Hz\n', min(rate));
fprintf('  Max:    %.4f bits/s/Hz\n', max(rate));
fprintf('\nSNR per User (mean across realizations):\n');
for k = 1:K
    fprintf('  User %d: %.2f ± %.2f dB\n', k, mean_snr(k), std_snr(k));
end
fprintf('\nIterations:\n');
fprintf('  Mean: %.1f\n', mean(iterations));
fprintf('  Max:  %d\n', max(iterations));
fprintf('========================================\n\n');

fprintf('Plotting completed successfully!\n');
