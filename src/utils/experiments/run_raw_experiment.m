function results = run_raw_experiment(config, project_root)
% RUN_RAW_EXPERIMENT Execute raw experiments and return results.
%
% Inputs:
%   config       - Struct with experiment settings
%   project_root - Path to project root (optional)
%
% Outputs:
%   results - Struct with raw experiment data

if nargin < 2
    project_root = '';
end

assert(isfield(config, 'num_trials'), 'Config missing field: num_trials');
assert(isfield(config, 'k_list'), 'Config missing field: k_list');
assert(isfield(config, 'scenarios'), 'Config missing field: scenarios');
assert(isfield(config, 'algorithms'), 'Config missing field: algorithms');

k_list = config.k_list(:)';
num_trials = config.num_trials;
scenarios = config.scenarios;
algorithms = config.algorithms;

if isfield(config, 'num_streams')
    num_streams = config.num_streams;
else
    num_streams = 1;
end

if isfield(config, 'noise_power')
    noise_power = config.noise_power;
else
    noise_power = 1.0;
end

if isfield(config, 'gamma_linear')
    gamma_linear = config.gamma_linear;
elseif isfield(config, 'gamma_db')
    gamma_linear = 10^(config.gamma_db / 10);
else
    error('Config must define gamma_linear or gamma_db');
end

channel_config = struct();
if isfield(config, 'channel')
    channel_config = config.channel;
end
channel_fn = resolve_channel(channel_config);

num_scenarios = numel(scenarios);
num_algs = numel(algorithms);
num_k = numel(k_list);

algo_names = cell(1, num_algs);
algo_fns = cell(1, num_algs);
method_fields = cell(1, num_algs);
for alg_idx = 1:num_algs
    algo_names{alg_idx} = algorithms(alg_idx).name;
    algo_fns{alg_idx} = resolve_algorithm(algorithms(alg_idx).name);
    method_fields{alg_idx} = matlab.lang.makeValidName(algo_names{alg_idx});
end
has_sdr = any(strcmpi(algo_names, 'sdr_beamformer'));

results = struct();
results.meta = struct();
results.meta.timestamp = datetime('now');
results.meta.seed = config.seed;
results.meta.project_root = project_root;

scenario_names = cell(1, num_scenarios);
for s = 1:num_scenarios
    if isfield(scenarios(s), 'name')
        scenario_names{s} = scenarios(s).name;
    else
        scenario_names{s} = sprintf('scenario_%d', s);
    end
end
results.meta.scenario_names = scenario_names;

for s = 1:num_scenarios
    scenario_name = scenario_names{s};
    scenario_field = matlab.lang.makeValidName(scenario_name);
    num_antennas = scenarios(s).num_antennas;

    results.raw.(scenario_field) = struct();
    results.raw.(scenario_field).name = scenario_name;
    results.raw.(scenario_field).num_antennas = num_antennas;
    results.raw.(scenario_field).k_list = k_list;
    results.raw.(scenario_field).algorithms = algo_names;
    results.raw.(scenario_field).methods = struct();
    for alg_idx = 1:num_algs
        method_field = method_fields{alg_idx};
        results.raw.(scenario_field).methods.(method_field) = struct();
        results.raw.(scenario_field).methods.(method_field).power = cell(1, num_k);
        results.raw.(scenario_field).methods.(method_field).min_snr = cell(1, num_k);
        results.raw.(scenario_field).methods.(method_field).snr = cell(1, num_k);
        results.raw.(scenario_field).methods.(method_field).feasible = cell(1, num_k);
        results.raw.(scenario_field).methods.(method_field).solve_time = cell(1, num_k);
        results.raw.(scenario_field).methods.(method_field).status = cell(1, num_k);
        results.raw.(scenario_field).methods.(method_field).cvx_optval = cell(1, num_k);
    end
    if has_sdr
        results.raw.(scenario_field).methods.BOUND = struct();
        results.raw.(scenario_field).methods.BOUND.cvx_optval = cell(1, num_k);
    end

    for k_idx = 1:num_k
        num_users = k_list(k_idx);

        if isscalar(gamma_linear)
            gamma = gamma_linear * ones(num_users, 1);
        else
            assert(numel(gamma_linear) == num_users, 'gamma_linear must be scalar or length K');
            gamma = gamma_linear(:);
        end

        if isscalar(noise_power)
            noise_k = noise_power * ones(num_users, 1);
        else
            assert(numel(noise_power) == num_users, 'noise_power must be scalar or length K');
            noise_k = noise_power(:);
        end

        for alg_idx = 1:num_algs
            method_field = method_fields{alg_idx};
            results.raw.(scenario_field).methods.(method_field).power{k_idx} = nan(1, num_trials);
            results.raw.(scenario_field).methods.(method_field).min_snr{k_idx} = nan(1, num_trials);
            results.raw.(scenario_field).methods.(method_field).snr{k_idx} = nan(num_users, num_trials);
            results.raw.(scenario_field).methods.(method_field).feasible{k_idx} = false(1, num_trials);
            results.raw.(scenario_field).methods.(method_field).solve_time{k_idx} = nan(1, num_trials);
            results.raw.(scenario_field).methods.(method_field).status{k_idx} = cell(1, num_trials);
            results.raw.(scenario_field).methods.(method_field).cvx_optval{k_idx} = nan(1, num_trials);
        end
        if has_sdr
            results.raw.(scenario_field).methods.BOUND.cvx_optval{k_idx} = nan(1, num_trials);
        end

        logf('Scenario %s (N=%d), K=%d', scenario_name, num_antennas, num_users);

        for trial = 1:num_trials
            H = channel_fn(num_antennas, num_users, num_streams);
            if ndims(H) == 3 && num_streams == 1
                H = squeeze(H);
            end

            for alg_idx = 1:num_algs
                method_field = method_fields{alg_idx};
                algo_params = struct();
                if isfield(algorithms(alg_idx), 'params') && ~isempty(algorithms(alg_idx).params)
                    algo_params = algorithms(alg_idx).params;
                end

                if isfield(algo_params, 'seed')
                    algo_params = rmfield(algo_params, 'seed');
                end

                algo_params.gamma = gamma;
                algo_params.noise_power = noise_k;

                try
                    [W, metrics] = algo_fns{alg_idx}(H, algo_params);

                    if isfield(metrics, 'snr') && ~isempty(metrics.snr)
                        snr = metrics.snr(:);
                    else
                        snr_config = struct('noise_power', noise_k);
                        snr = compute_snr(W, H, snr_config);
                    end

                    if isfield(metrics, 'final_power')
                        power = metrics.final_power;
                    else
                        power = norm(W, 'fro')^2;
                    end

                    min_snr = min(snr);

                    results.raw.(scenario_field).methods.(method_field).power{k_idx}(trial) = power;
                    results.raw.(scenario_field).methods.(method_field).min_snr{k_idx}(trial) = min_snr;
                    results.raw.(scenario_field).methods.(method_field).snr{k_idx}(:, trial) = snr;

                    if isfield(metrics, 'feasible')
                        feasible = metrics.feasible;
                    else
                        tol = 1e-6;
                        feasible = all(snr >= gamma * (1 - tol));
                    end
                    results.raw.(scenario_field).methods.(method_field).feasible{k_idx}(trial) = feasible;

                    if isfield(metrics, 'solve_time')
                        results.raw.(scenario_field).methods.(method_field).solve_time{k_idx}(trial) = metrics.solve_time;
                    end

                    if isfield(metrics, 'status_message')
                        results.raw.(scenario_field).methods.(method_field).status{k_idx}{trial} = metrics.status_message;
                    else
                        results.raw.(scenario_field).methods.(method_field).status{k_idx}{trial} = '';
                    end
                    if isfield(metrics, 'cvx_optval')
                        results.raw.(scenario_field).methods.(method_field).cvx_optval{k_idx}(trial) = metrics.cvx_optval;
                        if has_sdr && strcmpi(algo_names{alg_idx}, 'sdr_beamformer')
                            results.raw.(scenario_field).methods.BOUND.cvx_optval{k_idx}(trial) = metrics.cvx_optval;
                        end
                    end

                catch ME
                    results.raw.(scenario_field).methods.(method_field).status{k_idx}{trial} = ME.message;
                    logf('Algorithm %s failed (scenario=%s, K=%d, trial=%d): %s', ...
                        algo_names{alg_idx}, scenario_name, num_users, trial, ME.message);
                end
            end
        end
    end
end

end

function logf(fmt, varargin)
    if exist('log_message', 'file')
        log_message(fmt, varargin{:});
    else
        fprintf([fmt '\n'], varargin{:});
    end
end
