function [W, metrics] = zf_beamformer(H, config)
% ZF_BEAMFORMER Zero-Forcing baseline beamformer for multicast
%
% Implements the Zero-Forcing (ZF) approach as a baseline. The ZF beamformer
% nulls out inter-user interference by computing the pseudo-inverse, then
% scales to meet QoS requirements with minimum power.
%
% Inputs:
%   H      - [num_antennas x num_users] channel matrix
%            H(:,k) is the channel vector from transmitter to user k
%   config - struct with fields:
%            .gamma       : [num_users x 1] QoS/SNR targets (required)
%            .noise_power : [num_users x 1] noise power (required)
%
% Outputs:
%   W       - [num_antennas x 1] beamforming vector
%   metrics - struct with fields:
%             .final_power     : final transmit power
%             .snr             : [num_users x 1] SNR per user
%             .min_snr         : minimum SNR across users
%             .rate            : sum rate in bits/s/Hz
%             .feasible        : true if QoS constraints satisfied
%             .converged       : always true (closed-form solution)
%             .solve_time      : computation time in seconds
%             .status_message  : descriptive status string
%
% Method:
%   For multicast, all users receive the same signal. ZF computes a 
%   beamforming vector that points in the direction given by:
%   w = H * pinv(H' * H) * sqrt(P) * ones(K, 1)
%   Then scales to satisfy minimum SNR constraints.
%
% Reference:
%   Zero-forcing is a classical baseline in MIMO systems
%
% Note:
%   ZF requires num_antennas >= num_users for full rank

%% Configuration and parameter extraction

[num_antennas, num_users] = size(H);

% Extract QoS requirements (must be pre-generated vector)
if ~isfield(config, 'gamma')
    error('zf_beamformer:missingParameter', 'config.gamma is required');
end
gamma = config.gamma(:);
assert(length(gamma) == num_users, ...
    'gamma must be [num_users x 1] vector, got [%d x 1]', length(gamma));

% Extract noise power (must be pre-generated vector)
if ~isfield(config, 'noise_power')
    error('zf_beamformer:missingParameter', 'config.noise_power is required');
end
sigma_k_squared = config.noise_power(:);
assert(length(sigma_k_squared) == num_users, ...
    'noise_power must be [num_users x 1] vector, got [%d x 1]', length(sigma_k_squared));

%% Zero-Forcing beamformer computation

% For multicast beamforming, we want to transmit the same signal to all users
% ZF approach: Use pseudo-inverse to get initial direction
% w_zf = H * pinv(H' * H) * ones(K, 1)
%
% Alternative formulation for better numerical stability:
% w_zf = H * (H' * H)^(-1) * ones(K, 1)

% Check rank condition
if num_antennas < num_users
    % Return empty solution with NaN metrics
    W = zeros(num_antennas, 1);
    metrics = struct();
    metrics.final_power = NaN;
    metrics.snr = NaN(num_users, 1);
    metrics.min_snr = NaN;
    metrics.snr_db = NaN(num_users, 1);
    metrics.min_snr_db = NaN;
    metrics.power_db = NaN;
    metrics.rate = NaN;
    metrics.feasible = false;
    metrics.converged = false;
    metrics.solve_time = 0;
    metrics.cond_number = NaN;
    metrics.status_message = sprintf('ZF failed: underdetermined system (M=%d < K=%d)', num_antennas, num_users);
    return;
end

% Compute Gram matrix
G = H' * H;  % [num_users x num_users]

% Check condition number for numerical stability
cond_G = cond(G);
if cond_G > 1e10
    warning('zf_beamformer:illConditioned', ...
        'Channel Gram matrix is ill-conditioned (cond = %.2e). Results may be inaccurate.', cond_G);
end

% Compute ZF direction (for multicast: equal weight to all users)
try
    % Try direct inversion first
    w_zf = H / G * ones(num_users, 1);
catch
    % Fall back to pseudo-inverse if direct inversion fails
    w_zf = H * pinv(H' * H) * ones(num_users, 1);
end

% Normalize to unit norm (power will be scaled later)
w_zf = w_zf / norm(w_zf);

%% Scale to satisfy QoS constraints

% Compute received power at each user for unit-norm beamformer
recv_power_unit = abs(H' * w_zf).^2;  % [num_users x 1]

% Required power to satisfy each user's SNR constraint:
% SNR_k = |h_k' * w|^2 / sigma_k^2 >= gamma_k
% |h_k' * w|^2 >= gamma_k * sigma_k^2
%
% For scaled beamformer w = alpha * w_zf:
% alpha^2 * recv_power_unit(k) >= gamma_k * sigma_k^2
% alpha^2 >= gamma_k * sigma_k^2 / recv_power_unit(k)

alpha = compute_scaling_factor(w_zf, H, gamma, sigma_k_squared);

% Final beamformer
W = alpha * w_zf;


%% Compute metrics using centralized function
metrics = compute_beamformer_metrics(W, H, config);

% Add algorithm-specific fields
metrics.cond_number = cond_G;

% Add status message
if cond_G > 1e10
    metrics.status_message = sprintf('ZF (closed-form), WARNING: ill-conditioned (cond=%.2e)', cond_G);
elseif ~metrics.feasible
    metrics.status_message = 'WARNING: QoS constraints not satisfied at actual power (ZF)';
else
    metrics.status_message = 'ZF (closed-form solution)';
end

end
