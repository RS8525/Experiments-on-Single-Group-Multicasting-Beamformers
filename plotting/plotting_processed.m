function plotting_processed(processed_path, stat_name, metric_name, varargin)
% PLOTTING_PROCESSED Plot processed metrics vs K for multiple methods.
%
% Inputs:
%   processed_path - Path to processed_results.mat
%   stat_name      - 'mean', 'median', 'p25', 'p75', 'min', or 'max'
%   metric_name    - 'power', 'min_snr', 'snr_mean', 'cvx_optval', 'feasible', 'solve_time', 'power_db', 'min_snr_db'
%
% Name-value options:
%   'Scenario' - Scenario name or field (optional)
%   'Methods'  - Cell array of method names to plot (optional)
%   'SavePath' - Output file path for figure (optional)
%   'Title'    - Custom title (optional)

if nargin < 3
    error('Usage: plotting_processed(processed_path, stat_name, metric_name, ...)');
end

loaded = load(processed_path);
if ~isfield(loaded, 'processed')
    error('Expected variable "processed" in %s', processed_path);
end

processed = loaded.processed;

parser = inputParser;
addParameter(parser, 'Scenario', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Methods', {}, @(x) iscell(x) || isstring(x));
addParameter(parser, 'SavePath', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Title', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

scenario_name = char(parser.Results.Scenario);
method_filter = parser.Results.Methods;
save_path = char(parser.Results.SavePath);
title_text = char(parser.Results.Title);

stat_name = normalize_stat_name(stat_name);
metric_name = normalize_metric_name(metric_name);

scenario_fields = fieldnames(processed.scenarios);
if isempty(scenario_fields)
    error('No scenarios found in processed results');
end

if isempty(scenario_name)
    scenario_field = scenario_fields{1};
else
    scenario_field = matlab.lang.makeValidName(scenario_name);
    if ~isfield(processed.scenarios, scenario_field)
        error('Scenario not found: %s', scenario_name);
    end
end

scenario = processed.scenarios.(scenario_field);
k_list = scenario.k_list;

method_fields = fieldnames(scenario.methods);
if ~isempty(method_filter)
    method_fields = filter_methods(method_fields, method_filter);
end

figure('Position', [120, 120, 800, 500]);
hold on;
grid on;

line_styles = lines(numel(method_fields));
legend_entries = {};

for i = 1:numel(method_fields)
    method_field = method_fields{i};
    method = scenario.methods.(method_field);
    if ~isfield(method, 'metrics') || ~isfield(method.metrics, metric_name)
        fprintf('Skipping method %s (metric not found: %s)\n', method_field, metric_name);
        continue;
    end

    metric = method.metrics.(metric_name);
    if ~isfield(metric, stat_name)
        fprintf('Skipping method %s (stat not found: %s)\n', method_field, stat_name);
        continue;
    end

    y = metric.(stat_name);
    plot(k_list, y, 'LineWidth', 2, 'Color', line_styles(i, :));
    legend_entries{end+1} = format_method_name(method_field); %#ok<AGROW>
end

xlabel('Number of Users (K)');
ylabel(sprintf('%s (%s)', metric_label(metric_name), stat_name));

if isempty(title_text)
    title_text = sprintf('%s: %s', scenario.name, metric_label(metric_name));
end
title(title_text);

if ~isempty(legend_entries)
    legend(legend_entries, 'Location', 'best');
end

if ~isempty(save_path)
    saveas(gcf, save_path);
end
end

function stat_name = normalize_stat_name(stat_name)
    stat_name = char(lower(string(stat_name)));
    valid = {'mean', 'median', 'p25', 'p75', 'min', 'max'};
    if ~ismember(stat_name, valid)
        error('Invalid stat_name: %s', stat_name);
    end
end

function metric_name = normalize_metric_name(metric_name)
    metric_name = char(lower(string(metric_name)));
    if strcmp(metric_name, 'snr')
        metric_name = 'min_snr';
    end
    valid = {'power', 'min_snr', 'snr_mean', 'cvx_optval', 'feasible', 'solve_time', 'power_db', 'min_snr_db'};
    if ~ismember(metric_name, valid)
        error('Invalid metric_name: %s', metric_name);
    end
end

function label = metric_label(metric_name)
    switch metric_name
        case 'power'
            label = 'Transmit Power';
        case 'min_snr'
            label = 'Min SNR';
        case 'snr_mean'
            label = 'Mean SNR';
        case 'cvx_optval'
            label = 'SDR Bound (cvx optval)';
        case 'feasible'
            label = 'Feasibility Rate';
        case 'solve_time'
            label = 'Solve Time (s)';
        case 'power_db'
            label = metric_name;
        case 'min_snr_db'
            label = metric_name;
        otherwise
    end
end

function method_fields = filter_methods(method_fields, method_filter)
    if isstring(method_filter)
        method_filter = cellstr(method_filter);
    end
    normalized = cellfun(@(x) matlab.lang.makeValidName(char(x)), method_filter, 'UniformOutput', false);
    method_fields = intersect(method_fields, normalized, 'stable');
    if isempty(method_fields)
        error('No matching methods found');
    end
end

function name = format_method_name(method_field)
    if strcmpi(method_field, 'BOUND')
        name = '[BOUND]';
    else
        name = method_field;
    end
end
