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
%            .rc_c2_fn       : Reference constraint function handle (optional)
%
% Outputs:
%   W       - [num_antennas x 1] beamforming vector
%   metrics - struct with fields:
%             .num_iters      : number of iterations performed
%             .converged      : true if all constraints satisfied
%             .final_power    : final transmit power
%             .snr            : [num_users x 1] SNR per user
%             .min_snr        : minimum SNR across users
%             .rate           : sum rate in bits/s/Hz
%             .feasible       : true if QoS constraints satisfied
%             .h_orth         : orthogonalized channels
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

% Extract RC-C2 function (default: simple projection)
if isfield(config, 'rc_c2_fn')
    rc_c2_fn = config.rc_c2_fn;
else
    % Default orthogonalization function
    rc_c2_fn = @(h_k, h_l) h_k - (h_k' * h_l) / (h_l' * h_l) * h_l;
end

%% Main SBFC Algorithm

% Compute cost vector c = gamma .* sigma_k^2
c = gamma .* sigma_k_squared;

% Step 1: Save original channels and initialize orthogonalized channels
h_orth = H;

% Step 2: Select weakest user (most violated constraint)
[~, l1] = min(c .* vecnorm(h_orth, 2, 1).^2);

% Step 3: Set initial filter
w = h_orth(:, l1) / norm(h_orth(:, l1), 2);

% Step 4: Initialize active user set
U = setdiff(1:num_users, l1);

% Main iteration loop
converged = false;
snr_gamma_iter = ones(1, num_users);
l_prev = l1;

for iter = 2:num_antennas
    % Step 6: Find most violated constraint
    
    snr_gamma_iter(U) = abs(w' * H(:, U)).^2 ./ c(U);
    
    
    
    
    [min_snr_ratio, idx] = min(snr_gamma_iter(U));
    l_iter = U(idx);
    
    % Step 7: Exit if all constraints met
    if min_snr_ratio >= 1
        converged = true;
        break;
    end
    
    % Step 8: Orthogonalization
    h_norm_sq = vecnorm(h_orth(:, l_prev), 2)^2;
    h_orth(:, U) = h_orth(:, U) - (h_orth(:, l_prev)' * h_orth(:, U)) / h_norm_sq * h_orth(:, l_prev);
    
    % Step 9: Update beamformer
    alpha_iter = exp(-angle(w' * H(:, l_iter)) * 1i) * (sqrt(abs(w' * H(:, l_iter))^2 - c(l_iter)*snr_gamma_iter(l_iter)) - abs(w' * H(:, l_iter))) / vecnorm(h_orth(:, l_iter), 2);
    w = w + alpha_iter * h_orth(:, l_iter);
    
    % Step 10: Update active set
    U = setdiff(U, l_iter);
    
    % Exit if no more users in active set
    if isempty(U)
        converged = true;
        break;
    end
    l_prev = l_iter;
end




%% Compute final metrics
% Check feasability for num_users > num_antennas case and scale beamformer
alpha_final = compute_scaling_factor(w, H, gamma, sigma_k_squared);
W = alpha_final * w;
recv_power = abs(H' * W).^2;  

% SNR for each user
snr = recv_power ./ sigma_k_squared;

% Minimum SNR across users
min_snr = min(snr);

% Sum rate in bits/s/Hz
rate = sum(log2(1 + snr));

% Final transmit power
final_power = norm(W)^2;

% Check feasibility (with small tolerance)
tol_feas = 1e-6;
feasible = all(snr >= gamma * (1 - tol_feas));

% Populate metrics structure
metrics.num_iters = iter;
metrics.converged = converged;
metrics.final_power = final_power;
metrics.snr = snr;
metrics.min_snr = min_snr;
metrics.rate = rate;
metrics.feasible = feasible;
metrics.h_orth = h_orth;

if converged
    metrics.status_message = sprintf('Converged after %d iterations', iter);
else
    metrics.status_message = sprintf('Max iterations reached (%d)', max_iters);
end

end

