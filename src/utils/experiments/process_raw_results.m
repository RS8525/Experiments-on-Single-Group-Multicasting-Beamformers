function processed = process_raw_results(raw_results_path, output_path)
% PROCESS_RAW_RESULTS Compute summary statistics from raw experiment results.
%
% Inputs:
%   raw_results_path - Path to raw_results.mat
%   output_path      - (optional) Path to save processed_results.mat
%
% Outputs:
%   processed - Struct with processed statistics

if nargin < 1 || isempty(raw_results_path)
    error('raw_results_path is required');
end

if ~isfile(raw_results_path)
    error('Raw results file not found: %s', raw_results_path);
end

loaded = load(raw_results_path);
if isfield(loaded, 'results')
    results = loaded.results;
else
    error('Expected variable "results" in %s', raw_results_path);
end

processed = struct();
processed.meta = results.meta;
processed.meta.processed_timestamp = datetime('now');
processed.meta.raw_results_path = raw_results_path;
processed.config = results.meta.config_file;

scenario_fields = fieldnames(results.raw);
for s = 1:numel(scenario_fields)
    scenario_field = scenario_fields{s};
    scenario = results.raw.(scenario_field);

    processed.scenarios.(scenario_field) = struct();
    processed.scenarios.(scenario_field).name = scenario.name;
    processed.scenarios.(scenario_field).num_antennas = scenario.num_antennas;
    processed.scenarios.(scenario_field).k_list = scenario.k_list;
    processed.scenarios.(scenario_field).algorithms = scenario.algorithms;
    processed.scenarios.(scenario_field).methods = struct();

    method_fields = fieldnames(scenario.methods);
    num_k = numel(scenario.k_list);

    for m = 1:numel(method_fields)
        method_field = method_fields{m};
        method = scenario.methods.(method_field);

        metrics = struct();
        metrics.power = init_stats(num_k);
        metrics.min_snr = init_stats(num_k);
        metrics.snr_mean = init_stats(num_k);
        metrics.feasible = init_stats(num_k);
        metrics.solve_time = init_stats(num_k);
        metrics.power_db = init_stats(num_k);
        metrics.snr_db = init_stats(num_k);
        metrics.min_snr_db = init_stats(num_k);

        for k_idx = 1:num_k
            if isfield(method, 'power') && ~isempty(method.power)
                metrics.power = assign_stats(metrics.power, k_idx, method.power{k_idx});
            end
            if isfield(method, 'min_snr') && ~isempty(method.min_snr)
                metrics.min_snr = assign_stats(metrics.min_snr, k_idx, method.min_snr{k_idx});
            end
            if isfield(method, 'snr') && ~isempty(method.snr)
                snr_mat = method.snr{k_idx};
                if ~isempty(snr_mat)
                    snr_mean = mean(snr_mat, 1, 'omitnan');
                    metrics.snr_mean = assign_stats(metrics.snr_mean, k_idx, snr_mean);
                end
            end
            if isfield(method, 'feasible') && ~isempty(method.feasible)
                metrics.feasible = assign_stats(metrics.feasible, k_idx, double(method.feasible{k_idx}));
            end
            if isfield(method, 'solve_time') && ~isempty(method.solve_time)
                metrics.solve_time = assign_stats(metrics.solve_time, k_idx, method.solve_time{k_idx});
            end
            if isfield(method, 'power_db') && ~isempty(method.power_db)
                metrics.power_db = assign_stats(metrics.power_db, k_idx, method.power_db{k_idx});
            end
            if isfield(method, 'snr_db') && ~isempty(method.snr_db)
                snr_db_mat = method.snr_db{k_idx};
                if ~isempty(snr_db_mat)
                    snr_db_mean = mean(snr_db_mat, 1, 'omitnan');
                    metrics.snr_db = assign_stats(metrics.snr_db, k_idx, snr_db_mean);
                end
            end
            if isfield(method, 'min_snr_db') && ~isempty(method.min_snr_db)
                metrics.min_snr_db = assign_stats(metrics.min_snr_db, k_idx, method.min_snr_db{k_idx});
            end
            if isfield(method, 'exitflag') && ~isempty(method.exitflag)
                if ~isfield(metrics, 'exitflag_frequency')
                    metrics.exitflag_frequency = cell(1, num_k);
                end
                % Compute frequency distribution for exitflag values
                exitflag_values = method.exitflag{k_idx};
                exitflag_values = exitflag_values(~isnan(exitflag_values)); % Remove NaN
                if ~isempty(exitflag_values)
                    unique_flags = unique(exitflag_values);
                    frequencies = zeros(size(unique_flags));
                    for f = 1:length(unique_flags)
                        frequencies(f) = sum(exitflag_values == unique_flags(f)) / length(exitflag_values);
                    end
                    metrics.exitflag_frequency{k_idx} = struct('flags', unique_flags, 'frequencies', frequencies);
                else
                    metrics.exitflag_frequency{k_idx} = struct('flags', [], 'frequencies', []);
                end
            end
            if isfield(method, 'rank') && ~isempty(method.rank)
                if ~isfield(metrics, 'rank_frequency')
                    metrics.rank_frequency = cell(1, num_k);
                end
                % Compute frequency distribution for rank values
                rank_values = method.rank{k_idx};
                rank_values = rank_values(~isnan(rank_values)); % Remove NaN
                if ~isempty(rank_values)
                    unique_ranks = unique(rank_values);
                    frequencies = zeros(size(unique_ranks));
                    for r = 1:length(unique_ranks)
                        frequencies(r) = sum(rank_values == unique_ranks(r)) / length(rank_values);
                    end
                    metrics.rank_frequency{k_idx} = struct('ranks', unique_ranks, 'frequencies', frequencies);
                else
                    metrics.rank_frequency{k_idx} = struct('ranks', [], 'frequencies', []);
                end
            end
        end

        processed.scenarios.(scenario_field).methods.(method_field).metrics = metrics;
    end
end

if nargin < 2 || isempty(output_path)
    output_path = default_output_path(raw_results_path, processed);
end

output_dir = fileparts(output_path);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

save(output_path, 'processed');
end

function stats = init_stats(num_k)
    stats = struct();
    stats.mean = nan(1, num_k);
    stats.median = nan(1, num_k);
    stats.p25 = nan(1, num_k);
    stats.p75 = nan(1, num_k);
    stats.min = nan(1, num_k);
    stats.max = nan(1, num_k);
end

function stats = assign_stats(stats, k_idx, data)
    if isempty(data)
        return;
    end
    data = data(:);
    data = data(~isnan(data));
    if isempty(data)
        return;
    end
    stats.mean(k_idx) = mean(data);
    stats.median(k_idx) = median(data);
    q = prctile(data, [25 75]);
    stats.p25(k_idx) = q(1);
    stats.p75(k_idx) = q(2);
    stats.min(k_idx) = min(data);
    stats.max(k_idx) = max(data);
end

function output_path = default_output_path(raw_results_path, processed)
    raw_dir = fileparts(raw_results_path);
    results_root = fileparts(fileparts(raw_dir));

    experiment_name = '';
    if isfield(processed, 'meta') && isfield(processed.meta, 'experiment_name')
        experiment_name = processed.meta.experiment_name;
    end
    if isempty(experiment_name)
        [~, experiment_name] = fileparts(raw_dir);
    end

    output_dir = fullfile(results_root, 'processed', experiment_name);
    output_path = fullfile(output_dir, 'processed_results.mat');
end
