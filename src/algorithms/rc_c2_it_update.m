function [W, metrics] = rc_c2_it_update(H, config) 
% RC_C2_IT_UPDATE  Starts with RC-C2 and then performes the iterative SNR increasing update algorithm
% Implements the RC-C2 algorithm with an iterative update mechanism.
% Inputs:
%   H      - [num_antennas x num_users] channel matrix
%            H(:,k) is the channel vector from transmitter to user k
%   config - struct with fields:
%            .gamma       : [num_users x 1] QoS/SNR targets (required)
%            .noise_power : [num_users x 1] noise power (required) 
%            .max_iterations : maximum number of iterations (optional, default: 100)   
%            .alpha_max : maximum step size (optional, default: 0.2)
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
% Reference:
%   Paper: "Design of Single-Group Multicasting-Beamformers"


%% Configuration and defaults
[num_antennas, num_users] = size(H);

% Extract QoS requirements (must be pre-generated vector)
if ~isfield(config, 'gamma')
    error('rc_c2_it_update:missingParameter', 'config.gamma is required');
end
gamma = config.gamma(:);
assert(length(gamma) == num_users, ...
    'gamma must be [num_users x 1] vector, got [%d x 1]', length(gamma));

% Extract noise power (must be pre-generated vector)
if ~isfield(config, 'noise_power')
    error('rc_c2_it_update:missingParameter', 'config.noise_power is required');
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

% Extract alpha_max (maximum step size)
if isfield(config, 'alpha_max')
    alpha_max = config.alpha_max;
else
    alpha_max = 0.2;  % Default as per Algorithm 2
end

%% Initialize beamformer using RC-C2
% Use RC-C2 to get initial beamformer w^(0)
[w, ~] = rc_c2(H, config);

%% Initialize algorithm parameters
constraint_constant = gamma .* sigma_k_squared;  % Cost vector
alpha = alpha_max;  % Current step size
new_w = 1;          % Semaphore: 1 indicates w changed in previous iteration
dw = [];            % Update direction (saved between iterations)
recv_power_init = abs(H' * w).^2;
m_min = min(recv_power_init ./ constraint_constant); 

%% Main iterative update loop
for n = 1:max_iters
    
    % Lines 4-10: Compute update direction when w has changed
    if new_w == 1
        % Line 5: Find weakest constraint (search 'weakest' constr.)
        % Compute SNR for each user
        recv_power = abs(H' * w).^2;  % [K x 1]
        
        
        % Find user with minimum SNR/gamma ratio
        [m_min, ell] = min(recv_power ./ constraint_constant);
        
        % Line 7: Orthogonalize channel (Eq. 38)
        % h̃_ℓ(n) = h_ℓ(n) - (w^(n)^H * h_ℓ(n)) / P_tr * w^(n)
        h_ell = H(:, ell);
        P_tr = norm(w)^2;  % Current transmit power
        
        % Orthogonalization: remove component of w from h_ell
        h_ell_orth = h_ell - (w' * h_ell) / P_tr * w;
        
        % Line 8: Normalize to get update direction (Eq. 39)
        % w^(n)⊥ = h̃_ℓ(n) / ||h̃_ℓ(n)||_2
        norm_h_orth = norm(h_ell_orth);
        
        if norm_h_orth > eps
            w_perp = h_ell_orth / norm_h_orth;
        else
            % If orthogonalized channel is zero, use original direction
            w_perp = h_ell / norm(h_ell);
        end
        
        % Line 9: Save update direction
        dw = w_perp;
    end
    
    % Line 11: Compute magnitude of β (Eq. 32)
    % |β| = √(P_tr(2α - α²))
    beta_mag = sqrt(P_tr * (2*alpha - alpha^2));
    
    % Compute phase of β: angle(w^(n)^H * h_ℓ(n))
    % Note: We use the current w and the direction dw
    beta_phase = - angle(w' * h_ell);
    beta = beta_mag * exp(1j * beta_phase);
    
    % Line 12: Compute tentative update (Eq. 29)
    % w_temp = (1 - α)w^(n) + β·dw
    w_temp = (1 - alpha) * w + beta * dw;
    
    % Line 13: Evaluate if update improves minimum metric
    % Compute SNRs for tentative beamformer
    recv_power_temp = abs(H' * w_temp).^2;
    min_metric_temp = min(recv_power_temp ./ constraint_constant);
    
    % Lines 13-21: Accept or reject update
    if min_metric_temp > m_min
        % Lines 14-16: Accept update
        w = w_temp;           % Perform update
        new_w = 1;            % Indicate change of w
        
        % Line 16: Increase step-size if α ≤ α_max/2
        if alpha <= alpha_max / 2
            alpha = 2 * alpha;  % Increase step-size again
        end
        
    else
        % Lines 17-21: Reject update
        % w^(n+1) = w^(n)      (no beamformer change)
        % w stays the same
        new_w = 0;            % Indicate NO change of w
        
        % Line 20: Reduce step-size
        alpha = alpha / 2;
    end
end

%% Compute final metrics
% Final beamforming vector
W = w;

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
metrics.power_db = 10*log10(best_power);
metrics.snr_db = 10*log10(snr);
metrics.min_snr_db = 10*log10(min_snr);
metrics.num_iters = max_iters;
metrics.final_power = final_power;
metrics.snr = snr;
metrics.min_snr = min_snr;
metrics.rate = rate;
metrics.feasible = feasible;

if ~feasible
    metrics.status_message = [metrics.status_message, ' (WARNING: SNR constraints not fully satisfied)'];
end

end