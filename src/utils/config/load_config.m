function config = load_config(config_file)
% LOAD_CONFIG Load configuration from JSON file
%
% Inputs:
%   config_file - Path to JSON configuration file
%
% Outputs:
%   config - Struct containing configuration parameters
%
% Example:
%   config = load_config('../config/experiments/experiment_01.json');

% Check if file exists
if ~isfile(config_file)
    error('Configuration file not found: %s', config_file);
end

% Read JSON file
json_text = fileread(config_file);

% Parse JSON
try
    config = jsondecode(json_text);
catch ME
    error('Failed to parse JSON file: %s\n%s', config_file, ME.message);
end

fprintf('Configuration loaded from: %s\n', config_file);

end
