function gamma_k = generate_qos_parameters(num_users, config)
% GENERATE_QOS_PARAMETERS Generate QoS/SNR target parameters for all users
%
% This function generates the QoS (Quality of Service) target parameters
% (gamma_k) for all users based on the experiment configuration. It handles
% conversion from dB to linear scale and broadcasts scalar values to vectors.
%
% Inputs:
%   num_users - Number of users in the system
%   config    - Configuration struct with fields:
%               .gamma_linear : Scalar or [num_users x 1] linear QoS targets (optional)
%               .gamma_db     : Scalar or [num_users x 1] QoS targets in dB (optional)
%               .gamma_k      : Alias for gamma_linear (optional)
%               .gamma        : Alias for gamma_linear (optional)
%
% Output:
%   gamma_k   - [num_users x 1] vector of QoS/SNR targets in linear scale
%
% Notes:
%   - Priority order: gamma_linear > gamma_db > gamma_k > gamma
%   - If scalar value provided, it will be broadcast to all users
%   - At least one gamma field must be present in config
%   - gamma_db values are converted to linear: gamma_linear = 10^(gamma_db/10)
%
% Example:
%   config = struct('gamma_db', 0);  % 0 dB for all users
%   gamma_k = generate_qos_parameters(4, config);  % Returns [1; 1; 1; 1]

% Input validation
assert(isnumeric(num_users) && num_users > 0 && mod(num_users, 1) == 0, ...
    'num_users must be a positive integer');
assert(isstruct(config), 'config must be a struct');

% Extract gamma from config with priority order
if isfield(config, 'gamma_linear')
    gamma = config.gamma_linear;
elseif isfield(config, 'gamma_db')
    gamma = 10.^(config.gamma_db / 10);  % Convert dB to linear
elseif isfield(config, 'gamma_k')
    gamma = config.gamma_k;
elseif isfield(config, 'gamma')
    gamma = config.gamma;
else
    error('generate_qos_parameters:missingField', ...
        'Config must contain one of: gamma_linear, gamma_db, gamma_k, or gamma');
end

% Handle scalar or vector input
if isscalar(gamma)
    % Broadcast scalar to all users
    gamma_k = gamma * ones(num_users, 1);
else
    % Validate vector dimensions
    gamma = gamma(:);  % Ensure column vector
    if length(gamma) ~= num_users
        error('generate_qos_parameters:dimensionMismatch', ...
            'gamma must be scalar or have length num_users (%d), but got length %d', ...
            num_users, length(gamma));
    end
    gamma_k = gamma;
end

% Validate output
assert(all(gamma_k > 0), 'All QoS parameters must be positive');
assert(all(isfinite(gamma_k)), 'All QoS parameters must be finite');

end
