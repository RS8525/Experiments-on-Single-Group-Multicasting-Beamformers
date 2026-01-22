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
%            .max_iterations : Maximum iterations (optional, default = 100)
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
for iter = 2:max_iters
    % Step 6: Find most violated constraint
    snr_k = zeros(1, num_users);
    for k = U
        snr_k(k) = abs(w' * h_orth(:, k))^2;
    end
    
    snr_all = abs(w' * h_orth).^2;
    snr_ratio = snr_k(U) ./ snr_all(U);
    
    [~, idx] = max(snr_ratio);
    l_n = U(idx);
    
    % Step 7: Exit if all constraints met
    if snr_ratio(idx) >= 1
        converged = true;
        break;
    end
    
    % Step 8: Orthogonalization
    for k = U
        h_orth(:, k) = rc_c2_fn(h_orth(:, k), h_orth(:, l_n));
    end
    
    % Step 9: Update filter
    alpha_n = compute_alpha(snr_k(l_n), snr_all(l_n), w, h_orth(:, l_n));
    w = w + alpha_n * (h_orth(:, l_n) / norm(h_orth(:, l_n), 2));
    
    % Step 10: Update active set
    U = setdiff(U, l_n);
    
    % Exit if no more users in active set
    if isempty(U)
        converged = true;
        break;
    end
end

% Output beamforming vector
W = w;

%% Compute final metrics
% Received power at each user
recv_power = abs(H' * W).^2;  % [num_users x 1]

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

%% Helper function
function alpha = compute_alpha(snr_k, snr_l, w, h_l)
% COMPUTE_ALPHA Compute step size for filter update
%
% Inputs:
%   snr_k : SNR at user k
%   snr_l : Total SNR at user l
%   w     : Current beamforming vector
%   h_l   : Channel vector for user l
%
% Output:
%   alpha : Step size for beamforming update

h_l_norm = norm(h_l, 2);
if h_l_norm < eps
    alpha = 0;
else
    alpha = (sqrt(snr_k) - sqrt(snr_l)) / h_l_norm;
end

end