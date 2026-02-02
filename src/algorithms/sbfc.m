function [W, metrics] = sbfc(H, config)
% SBFC Successive Beamforming Filter Computation
%
% Implements the SBFC algorithm for multicast beamforming using successive
% constraint satisfaction and channel orthogonalization.
%
% Inputs:
%   H      - [num_antennas x num_users] channel matrix
%            H(:,k) is the channel vector from transmitter to user k
%   config - struct with fields:
%            .gamma          : [num_users x 1] QoS/SNR targets (required)
%            .noise_power    : [num_users x 1] noise power (required)
%
% Outputs:
%   W       - [num_antennas x 1] beamforming vector
%   metrics - struct with fields:
%             .num_iters      : number of iterations performed
%             .final_power    : final transmit power
%             .snr            : [num_users x 1] SNR per user
%             .min_snr        : minimum SNR across users
%             .rate           : sum rate in bits/s/Hz
%             .feasible       : true if QoS constraints satisfied
%             .status_message : descriptive status string

%% Configuration and defaults
[num_antennas, num_users] = size(H);

% Extract QoS requirements (must be pre-generated vector)
if ~isfield(config, 'gamma')
    error('sbfc:missingParameter', 'config.gamma is required');
end
gamma = config.gamma(:);
assert(length(gamma) == num_users, ...
    'gamma must be [num_users x 1] vector, got [%d x 1]', length(gamma));

% Extract noise power (must be pre-generated vector)
if ~isfield(config, 'noise_power')
    error('sbfc:missingParameter', 'config.noise_power is required');
end
sigma_k_squared = config.noise_power(:);
assert(length(sigma_k_squared) == num_users, ...
    'noise_power must be [num_users x 1] vector, got [%d x 1]', length(sigma_k_squared));

% Extract maximum iterations
if isfield(config, 'max_iterations')
    max_iters = config.max_iterations;
else
    max_iters = 100;
end


%% Main SBFC Algorithm

% Compute cost vector c = gamma .* sigma_k^2
c = gamma .* sigma_k_squared;

% Step 1: Save original channels and initialize orthogonalized channels
h_orth = H;

% Step 2: Select weakest user (most violated constraint)
[~, l1] = min((vecnorm(h_orth).^2)' ./ c);

% Step 3: Set initial filter
w = h_orth(:, l1) / norm(h_orth(:, l1), 2)^2 * sqrt(c(l1));

% Step 4: Initialize active user set
U = setdiff(1:num_users, l1);

% Main iteration loop
converged = false;
snr_gamma_iter = ones(1, num_users);
l_prev = l1;

for iter = 2:num_antennas
    % Step 6: Find most violated constraint
    
    % Compute SNR/gamma ratio for active users
    recv_power_U = abs(w' * H(:, U)).^2;
    snr_gamma_temp = recv_power_U' ./ c(U);
    
    [min_snr_ratio, idx] = min(snr_gamma_temp);
    snr_gamma_iter(U) = snr_gamma_temp(:);
    l_iter = U(idx);
    
    % Step 7: Exit if all constraints met
    if min_snr_ratio >= 1
        converged = true;
        break;
    end
    
    % Step 8: Orthogonalization
    h_norm_sq = vecnorm(h_orth(:, l_prev), 2)^2;
    proj_coeffs = (h_orth(:, l_prev)' * h_orth(:, U)) / h_norm_sq;
    h_orth(:, U) = h_orth(:, U) - h_orth(:, l_prev) * proj_coeffs;
    
    % Step 9: Update beamformer
    alpha_iter_angle = exp(-angle(w' * H(:, l_iter)) * 1i);
    alpha_iter_abs = (sqrt(c(l_iter) ) - abs(w' * H(:, l_iter))) / vecnorm(h_orth(:, l_iter), 2)^2;
    %Warning sigma_k_squared used in the formula under the square root but in the paper it is sigma_\mu_squared
    alpha_iter = alpha_iter_angle * alpha_iter_abs;
    w = w + alpha_iter * h_orth(:, l_iter);
    
    % Step 10: Update active set
    U = setdiff(U, l_iter);
    
    % Exit if no more users in active set
    if isempty(U)
        break;
    end
    l_prev = l_iter;
end




%% Compute final metrics
% Check feasability for num_users > num_antennas case and scale beamformer
alpha_final = compute_scaling_factor(w, H, gamma, sigma_k_squared);
W = alpha_final * w;

% Use centralized metrics computation for consistency
metrics = compute_beamformer_metrics(W, H, config);

% Add algorithm-specific metrics
metrics.num_iters = iter;

% Add algorithm-specific status message if infeasible
if ~metrics.feasible
    metrics.status_message = 'WARNING: QoS constraints not satisfied at actual power';
end

end

