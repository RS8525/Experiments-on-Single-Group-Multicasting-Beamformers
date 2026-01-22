function [W, metrics] = rc_c2(H, config)
% RC_C2 Reduced Complexity Combine-2 algorithm for multicast beamforming
%
% Implements the RC-C2 algorithm based on the heuristic that the user with
% the smallest metric c_k^(-1)||h_k||_2^2 is likely to be part of the optimal
% two-user set. This reduces complexity from K(K-1)/2 to K-1 pairs.
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
%             .num_pairs       : number of user pairs considered (K-1)
%             .best_pair       : indices [i, j] of the best user pair
%             .selected_user   : index of the heuristically selected user
%             .final_power     : final transmit power
%             .snr             : [num_users x 1] SNR per user
%             .min_snr         : minimum SNR across users
%             .rate            : sum rate in bits/s/Hz
%             .feasible        : true if QoS constraints satisfied
%             .status_message  : descriptive status string
%
% Reference:
%   Paper: "Design of Single-Group Multicasting-Beamformers"
%   Algorithm: Reduced Complexity Combine-2 (RC-C2)
%   Complexity: 4NK + 21K - 2N + O(1) FLOPs
    
%% Configuration and defaults
[num_antennas, num_users] = size(H);

% Extract QoS requirements (must be pre-generated vector)
if ~isfield(config, 'gamma')
    error('rc_c2:missingParameter', 'config.gamma is required');
end
gamma = config.gamma(:);
assert(length(gamma) == num_users, ...
    'gamma must be [num_users x 1] vector, got [%d x 1]', length(gamma));

% Extract noise power (must be pre-generated vector)
if ~isfield(config, 'noise_power')
    error('rc_c2:missingParameter', 'config.noise_power is required');
end
sigma_k_squared = config.noise_power(:);
assert(length(sigma_k_squared) == num_users, ...
    'noise_power must be [num_users x 1] vector, got [%d x 1]', length(sigma_k_squared));

%% Main algorithm: RC-C2 heuristic-based approach
% Heuristic: The user with the smallest metric c_k^(-1) * ||h_k||_2^2
% is likely to be part of the two-user set {i,j}
% This reduces complexity from K(K-1)/2 to K-1 pairs

% Compute constants c_k = gamma_k * sigma_k^2
c_k = gamma .* sigma_k_squared;

% Compute metric c_k^(-1) * ||h_k||_2^2 for each user
metric_k = zeros(num_users, 1);
for k = 1:num_users
    metric_k(k) = norm(H(:, k))^2 / c_k(k);
end

% Find the user with the smallest metric
[~, k_star] = min(metric_k);

% Loop over K-1 user pairs: pair user k_star with each other user
num_pairs = num_users - 1;
best_power = inf;
best_w = [];
best_pair = [0, 0];

for k = 1:num_users
    if k == k_star
        continue;  % Skip pairing k_star with itself
    end
    
    % Extract channel vectors for users k_star and k
    h_i = H(:, k_star);
    h_j = H(:, k);
    
    % Compute the lowest norm precoding vector w^{k_star,k}
    w_ij = compute_lowest_norm_precoder(h_i, h_j, gamma(k_star), gamma(k), ...
                                       sigma_k_squared(k_star), sigma_k_squared(k));
    
    % Compute scaling factor alpha^{k_star,k} to satisfy ALL SNR constraints
    alpha_ij = compute_scaling_factor(w_ij, H, gamma, sigma_k_squared);
    
    % Scale the precoder if necessary
    if alpha_ij > 1
        w_scaled = alpha_ij * w_ij;
    else
        w_scaled = w_ij;
    end
    
    % Compute transmit power for this candidate
    power_ij = norm(w_scaled)^2;
    
    % Keep track of the best (minimum power) solution
    if power_ij < best_power
        best_power = power_ij;
        best_w = w_scaled;
        best_pair = [k_star, k];
    end
end

% Select the final beamforming vector
W = best_w;

%% Compute final metrics
% Received power at each user
recv_power = abs(H' * W).^2;  % [num_users x 1]

% SNR for each user
snr = recv_power ./ sigma_k_squared;

% Minimum SNR across users
min_snr = min(snr);

% Sum rate in bits/s/Hz
rate = sum(log2(1 + snr));

% Check feasibility (with small tolerance)
feasible = all(snr >= gamma * (1 - 1e-8));

% Populate metrics structure
metrics.num_pairs = num_pairs;
metrics.best_pair = best_pair;
metrics.selected_user = k_star;
metrics.final_power = best_power;
metrics.snr = snr;
metrics.min_snr = min_snr;
metrics.rate = rate;
metrics.feasible = feasible;

if feasible
    metrics.status_message = sprintf('Solution found (pair: %d, %d; selected user: %d)', ...
                                     best_pair(1), best_pair(2), k_star);
else
    metrics.status_message = sprintf('Solution found but infeasible (pair: %d, %d; selected user: %d)', ...
                                     best_pair(1), best_pair(2), k_star);
end

end
