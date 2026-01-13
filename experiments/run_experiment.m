function results = run_experiment(config_path)
% RUN_EXPERIMENT Run raw experiments defined by a JSON config file.
%
% Inputs:
%   config_path - Path to JSON configuration file
%
% Outputs:
%   results - Struct with raw experiment data

if nargin < 1 || isempty(config_path)
    error('config_path is required');
end

% Add paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..');
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'config')));

% Resolve and load configuration
config_file = resolve_config_path(config_path, project_root);
config = load_config(config_file);

% Merge defaults if available
defaults_file = fullfile(project_root, 'config', 'defaults', 'default_params.json');
if isfile(defaults_file)
    defaults = load_config(defaults_file);
    config = merge_configs(defaults, config);
end

if ~isfield(config, 'experiment_name') || isempty(config.experiment_name)
    error('Config missing field: experiment_name');
end
if ~isfield(config, 'seed')
    error('Config missing field: seed');
end

% Set global seed once for reproducibility
set_global_seed(config.seed);

% Run experiment
results = run_raw_experiment(config, project_root);

% Add metadata
if ~isfield(results, 'meta') || isempty(results.meta)
    results.meta = struct();
end
results.meta.experiment_name = config.experiment_name;
results.meta.config_file = config_file;
results.meta.project_root = project_root;

% Save results (overwrite)
output_dir = fullfile(project_root, 'results', 'raw', config.experiment_name);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

output_file = fullfile(output_dir, 'raw_results.mat');
save(output_file, 'results');

logf('Results saved to: %s', output_file);

end

function config_file = resolve_config_path(config_path, project_root)
    if isfile(config_path)
        config_file = config_path;
        return;
    end

    candidate = fullfile(project_root, config_path);
    if isfile(candidate)
        config_file = candidate;
        return;
    end

    error('Config file not found: %s', config_path);
end

function logf(fmt, varargin)
    if exist('log_message', 'file')
        log_message(fmt, varargin{:});
    else
        fprintf([fmt '\n'], varargin{:});
    end
end
