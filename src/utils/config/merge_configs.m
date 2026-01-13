function merged = merge_configs(default_config, user_config)
% MERGE_CONFIGS Merge user config with default config
%
% Inputs:
%   default_config - Struct with default parameter values
%   user_config    - Struct with user-specified parameters
%
% Outputs:
%   merged - Struct with merged configuration (user values override defaults)
%
% Note: Only top-level fields are merged. Nested structs in user_config
%       completely replace the corresponding nested struct in default_config.
%
% Example:
%   defaults = load_config('../config/defaults/default_params.json');
%   user = load_config('../config/experiments/experiment_01.json');
%   config = merge_configs(defaults, user);

% Start with default config
merged = default_config;

% Override with user config fields
if ~isempty(user_config)
    user_fields = fieldnames(user_config);
    for i = 1:length(user_fields)
        field = user_fields{i};
        merged.(field) = user_config.(field);
    end
end

end
