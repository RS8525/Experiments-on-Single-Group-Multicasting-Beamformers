function [W, metrics] = ff_c2(H, config)
% FF_C2 Full Featured Combine-2 algorithm for multicast beamforming
%
% Implements the FF-C2 algorithm that selects two users out of K, computes
% the lowest norm precoding vector meeting their SNR constraints, scales it
% to satisfy all users' constraints, and chooses the solution with minimum
% transmit power.
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
%             .status_message  : descriptive status string
%
% Reference:
%   Paper: "Design of Single-Group Multicasting-Beamformers"
%   Algorithm: Full Featured Combine-2 (FF-C2)
%   Complexity: K^2(N+10) + KN + O(1) FLOPs

%% Configuration and defaults
[num_antennas, num_users] = size(H);

% Extract QoS requirements (must be pre-generated vector)
if ~isfield(config, 'gamma')
    error('ff_c2:missingParameter', 'config.gamma is required');
end
gamma = config.gamma(:);
assert(length(gamma) == num_users, ...
    'gamma must be [num_users x 1] vector, got [%d x 1]', length(gamma));

% Extract noise power (must be pre-generated vector)
if ~isfield(config, 'noise_power')
    error('ff_c2:missingParameter', 'config.noise_power is required');
end
sigma_k_squared = config.noise_power(:);
assert(length(sigma_k_squared) == num_users, ...
    'noise_power must be [num_users x 1] vector, got [%d x 1]', length(sigma_k_squared));

%% Main algorithm: Loop over all K(K-1)/2 user pairs
num_pairs = num_users * (num_users - 1) / 2;
best_power = inf;
best_w = [];
best_pair = [0, 0];

pair_idx = 0;
for i = 1:(num_users-1)
    for j = (i+1):num_users
        pair_idx = pair_idx + 1;
        
        % Extract channel vectors for users i and j
        h_i = H(:, i);
        h_j = H(:, j);
        
        % Compute the lowest norm precoding vector w^{i,j} for users i and j
        w_ij = compute_lowest_norm_precoder(h_i, h_j, gamma(i), gamma(j), ...
                                           sigma_k_squared(i), sigma_k_squared(j));
        
        % Compute scaling factor alpha^{i,j} to satisfy ALL SNR constraints
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
            best_pair = [i, j];
        end
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
metrics.power_db = 10*log10(best_power);
metrics.snr_db = 10*log10(snr);
metrics.min_snr_db = 10*log10(min_snr);
metrics.final_power = best_power;
metrics.snr = snr;
metrics.min_snr = min_snr;
metrics.rate = rate;
metrics.feasible = feasible;

if ~feasible
    metrics.status_message = [metrics.status_message, ' (WARNING: SNR constraints not fully satisfied)'];
end

end
