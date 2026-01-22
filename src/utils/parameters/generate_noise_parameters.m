function sigma_k_squared = generate_noise_parameters(num_users, config)
% GENERATE_NOISE_PARAMETERS Generate noise power parameters for all users
%
% This function generates the noise power parameters (sigma_k^2) for all
% users based on the experiment configuration. It handles broadcasting
% scalar values to vectors.
%
% Inputs:
%   num_users - Number of users in the system
%   config    - Configuration struct with fields:
%               .noise_power       : Scalar or [num_users x 1] noise power (optional, default=1.0)
%               .sigma_k_squared   : Alias for noise_power (optional)
%               .sigma_squared     : Alias for noise_power (optional)
%
% Output:
%   sigma_k_squared - [num_users x 1] vector of noise power values
%
% Notes:
%   - Priority order: noise_power > sigma_k_squared > sigma_squared > default (1.0)
%   - If scalar value provided, it will be broadcast to all users
%   - Default noise power is 1.0 if no field is present
%   - All noise power values must be positive and finite
%
% Example:
%   config = struct('noise_power', 1.0);
%   sigma_k_squared = generate_noise_parameters(4, config);  % Returns [1; 1; 1; 1]

% Input validation
assert(isnumeric(num_users) && num_users > 0 && mod(num_users, 1) == 0, ...
    'num_users must be a positive integer');
assert(isstruct(config), 'config must be a struct');

% Extract noise power from config with priority order
if isfield(config, 'noise_power')
    noise_power = config.noise_power;
elseif isfield(config, 'sigma_k_squared')
    noise_power = config.sigma_k_squared;
elseif isfield(config, 'sigma_squared')
    noise_power = config.sigma_squared;
else
    % Default noise power
    noise_power = 1.0;
end

% Handle scalar or vector input
if isscalar(noise_power)
    % Broadcast scalar to all users
    sigma_k_squared = noise_power * ones(num_users, 1);
else
    % Validate vector dimensions
    noise_power = noise_power(:);  % Ensure column vector
    if length(noise_power) ~= num_users
        error('generate_noise_parameters:dimensionMismatch', ...
            'noise_power must be scalar or have length num_users (%d), but got length %d', ...
            num_users, length(noise_power));
    end
    sigma_k_squared = noise_power;
end

% Validate output
assert(all(sigma_k_squared > 0), 'All noise power values must be positive');
assert(all(isfinite(sigma_k_squared)), 'All noise power values must be finite');

end
