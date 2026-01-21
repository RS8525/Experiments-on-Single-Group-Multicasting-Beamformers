% plot_example_results.m
% Example plotting script for multicast beamforming results (new schema)
%
% This script demonstrates:
% - Loading saved results from results/raw/<experiment_name>/raw_results.mat
% - Selecting scenario and method
% - Generating performance figures
% - Saving figures to results/figures/

% Add paths
addpath(genpath('../src'));

% Clear workspace
clear; close all; clc;

%% Load results
fprintf('Loading results...\n');

% Specify result file (modify this path to your actual result file)
result_file = '../results/raw/fig4/raw_results.mat';

% If the file doesn't exist, try to find the most recent one
if ~isfile(result_file)
    result_dir = '../results/raw';
    files = dir(fullfile(result_dir, '*', 'raw_results.mat'));

    if isempty(files)
        error('No raw_results.mat files found in %s', result_dir);
    end

    % Sort by date and get most recent
    [~, idx] = sort([files.datenum], 'descend');
    result_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    fprintf('Using most recent result file: %s\n', result_file);
end

load(result_file, 'results');

%% Select scenario and method
scenario_name = '';
method_name = '';
k_idx = 1;

scenario_fields = fieldnames(results.raw);
if isempty(scenario_fields)
    error('No scenarios found in results.raw');
end

if isempty(scenario_name)
    scenario_field = scenario_fields{1};
else
    scenario_field = matlab.lang.makeValidName(scenario_name);
    if ~isfield(results.raw, scenario_field)
        error('Scenario not found: %s', scenario_name);
    end
end

scenario = results.raw.(scenario_field);
method_name = select_method(scenario.methods, method_name);
method = scenario.methods.(method_name);

if k_idx < 1 || k_idx > numel(scenario.k_list)
    error('k_idx out of range. Available K count: %d', numel(scenario.k_list));
end

snr = method.snr{k_idx};
if isempty(snr)
    error('Selected method has no SNR data: %s', method_name);
end

snr_db = 10 * log10(snr);
rate = sum(log2(1 + snr), 1);
power = method.power{k_idx};
min_snr = method.min_snr{k_idx};

K = size(snr_db, 1);
num_trials = size(snr_db, 2);

%% Create output directory
output_dir = '../results/figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Figure 1: Distributions
figure('Position', [100, 100, 900, 600]);

subplot(2, 2, 1);
histogram(rate, 20, 'Normalization', 'probability');
xlabel('Sum Rate (bits/s/Hz)');
ylabel('Probability');
title('Sum Rate Distribution');
grid on;

subplot(2, 2, 2);
plot(rate, 'b.-', 'LineWidth', 1.5);
xlabel('Trial');
ylabel('Sum Rate (bits/s/Hz)');
title('Sum Rate vs Trial');
grid on;

subplot(2, 2, 3);
histogram(power, 20, 'Normalization', 'probability');
xlabel('Transmit Power');
ylabel('Probability');
title('Power Distribution');
grid on;

subplot(2, 2, 4);
histogram(min_snr, 20, 'Normalization', 'probability');
xlabel('Min SNR');
ylabel('Probability');
title('Min SNR Distribution');
grid on;

sgtitle(sprintf('%s - %s (N=%d, K=%d)', ...
    scenario.name, method_name, scenario.num_antennas, scenario.k_list(k_idx)));

fig_file = fullfile(output_dir, 'example_summary.png');
saveas(gcf, fig_file);
fprintf('Figure saved: %s\n', fig_file);

%% Figure 2: SNR per User
figure('Position', [150, 150, 800, 500]);

mean_snr = mean(snr_db, 2);
std_snr = std(snr_db, 0, 2);

errorbar(1:K, mean_snr, std_snr, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('User Index');
ylabel('SNR (dB)');
title('Mean SNR per User with Standard Deviation');
grid on;
ylim([min(mean_snr - std_snr) - 2, max(mean_snr + std_snr) + 2]);

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

median_rate = median(rate);
hold on;
plot([median_rate, median_rate], [0, 1], 'r--', 'LineWidth', 1.5);
legend('CDF', sprintf('Median (%.3f)', median_rate), 'Location', 'southeast');

fig_file = fullfile(output_dir, 'example_rate_cdf.png');
saveas(gcf, fig_file);
fprintf('Figure saved: %s\n', fig_file);

%% Display statistics
fprintf('\n========== Statistics Summary ==========\n');
fprintf('Scenario: %s\n', scenario.name);
fprintf('Method:   %s\n', method_name);
fprintf('K:        %d\n', scenario.k_list(k_idx));
fprintf('Trials:   %d\n', num_trials);
fprintf('\nSum Rate:\n');
fprintf('  Mean:   %.4f bits/s/Hz\n', mean(rate));
fprintf('  Median: %.4f bits/s/Hz\n', median(rate));
fprintf('  Std:    %.4f bits/s/Hz\n', std(rate));
fprintf('  Min:    %.4f bits/s/Hz\n', min(rate));
fprintf('  Max:    %.4f bits/s/Hz\n', max(rate));
fprintf('\nPower:\n');
fprintf('  Mean: %.4f\n', mean(power));
fprintf('  Std:  %.4f\n', std(power));
fprintf('  Min:  %.4f\n', min(power));
fprintf('  Max:  %.4f\n', max(power));
fprintf('\nSNR per User (mean across trials):\n');
for k = 1:K
    fprintf('  User %d: %.2f +- %.2f dB\n', k, mean_snr(k), std_snr(k));
end
fprintf('========================================\n\n');

fprintf('Plotting completed successfully!\n');

function method_name = select_method(methods, requested)
    if nargin >= 2 && ~isempty(requested)
        requested_field = matlab.lang.makeValidName(requested);
        if isfield(methods, requested_field)
            method_name = requested_field;
            return;
        end
        error('Method not found: %s', requested);
    end

    method_names = fieldnames(methods);
    for i = 1:numel(method_names)
        candidate = method_names{i};
        if isfield(methods.(candidate), 'snr') && ~isempty(methods.(candidate).snr)
            method_name = candidate;
            return;
        end
    end

    method_name = method_names{1};
end
