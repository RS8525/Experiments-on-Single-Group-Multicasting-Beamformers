function [W, metrics] = sdr_beamformer(H, config)
% SDR_BEAMFORMER Semi-definite relaxation beamformer without the rank = 1 constraint, after randomization is Applied
%
% Solves the multicast beamforming problem using semi-definite relaxation
% (rank relaxed) followed by Gaussian randomization to extract rank-1 solution.
%
% Inputs:
%   H      - [num_antennas x num_users] channel matrix
%            H(:,k) is the channel vector from transmitter to user k
%   config - struct with fields:
%            .gamma            : scalar or [num_users x 1] QoS requirements (optional, default = ones)
%            .noise_power        : scalar or [num_users x 1] noise power (optional, default = ones)
%            .num_randomizations : number of Gaussian randomizations (optional, default = 200)
%            .cvx_quiet          : suppress CVX output (optional, default = true)
%            .seed               : random seed for reproducibility (optional)
%
% Outputs:
%   W       - [num_antennas x 1] beamforming vector
%   metrics - struct with fields:
%             .cvx_status         : CVX solver status
%             .cvx_optval         : optimal value from CVX (trace of X)
%             .solve_time         : solver time in seconds
%             .best_power         : transmit power of best randomized solution
%             .final_power        : final transmit power (same as best_power)
%             .snr                : [num_users x 1] SNR per user
%             .min_snr            : minimum SNR across users
%             .rate               : sum rate in bits/s/Hz
%             .feasible           : true if QoS constraints satisfied
%             .converged          : true if CVX solved successfully
%             .status_message     : descriptive status string
%
% Requirements:
%   CVX toolbox must be installed
%
% Method:
%   1. Solve SDR: minimize trace(X) s.t. trace(Q_k*X) >= gamma*noise_k, X >= 0
%   2. Apply Gaussian randomization to extract rank-1 beamformer
%   3. Scale randomized solutions to satisfy QoS constraints with minimum power
%
% Reference:
%   Paper: "Design of Single-Group Multicasting-Beamformers"
%   Paper: "Transmit beamforming for physical-layer multicasting


%% Configuration and defaults

[num_antennas, num_users] = size(H);

% Set random seed if provided
if isfield(config, 'seed')
    % Try to use the project's set_global_seed utility
    if exist('set_global_seed', 'file')
        set_global_seed(config.seed);
    else
        rng(config.seed, 'twister');
    end
end

% Extract parameters with defaults
% Note: As of the refactoring, gamma and noise_power should be pre-generated
% vectors. However, we maintain backward compatibility by allowing defaults.
if isfield(config, 'gamma')
    gamma = config.gamma(:);
    assert(length(gamma) == num_users, ...
        'gamma must be [num_users x 1] vector, got [%d x 1]', length(gamma));
else
    warning('sdr_beamformer:missingGamma', ...
        'config.gamma not provided, using default value of 1.0 for all users');
    gamma = ones(num_users, 1);
end

if isfield(config, 'noise_power')
    noise_power = config.noise_power(:);
    assert(length(noise_power) == num_users, ...
        'noise_power must be [num_users x 1] vector, got [%d x 1]', length(noise_power));
else
    warning('sdr_beamformer:missingNoisePower', ...
        'config.noise_power not provided, using default value of 1.0 for all users');
    noise_power = ones(num_users, 1);
end

if isfield(config, 'num_randomizations')    
    num_randomizations = config.num_randomizations;
    assert(num_randomizations > 0, 'num_randomizations must be positive');
else
    num_randomizations = 100;  % Default
end


if isfield(config, 'cvx_quiet')
    cvx_quiet_mode = config.cvx_quiet;
else
    cvx_quiet_mode = true;  % Default to quiet
end

%% Precompute constants outside CVX
% Q{k} = h_k * h_k^H for each user
Q = cell(num_users, 1);
for k = 1:num_users
    h_k = H(:, k);
    Q{k} = h_k * h_k';  % [num_antennas x num_antennas]
end

% Right-hand side of QoS constraints: noise_k * gamma
rhs = noise_power .* gamma;  % [num_users x 1]

%% Solve SDR via CVX using SeDuMi

tic;
cvx_begin sdp
    if cvx_quiet_mode
        cvx_quiet(true);
    end
    
    variable X(num_antennas, num_antennas) hermitian semidefinite
    
    minimize( trace(X) )
    
    subject to
        for k = 1:num_users
            trace(Q{k} * X) >= rhs(k);  %#ok<VUNUS>
        end
cvx_end
solve_time = toc;

% Store CVX metrics
metrics.cvx_status = cvx_status;
metrics.solve_time = solve_time;

% Check if CVX solved successfully
if ~strcmp(cvx_status, 'Solved') && ~strcmp(cvx_status, 'Inaccurate/Solved')
    % CVX failed - return empty solution
    W = zeros(num_antennas, 1);
    metrics.cvx_optval = inf;
    metrics.num_randomizations = 0;
    metrics.best_power = inf;
    metrics.final_power = inf;
    metrics.snr = zeros(num_users, 1);
    metrics.min_snr = 0;
    metrics.rate = 0;
    metrics.feasible = false;
    metrics.converged = false;
    metrics.status_message = sprintf('CVX failed: %s', cvx_status);
    return;
end

X_opt = X;



%% Gaussian randomization to extract rank-1 beamformer
if ~(X_opt == X_opt')  % Ensure Hermitian
    printf('Error: X_opt is not Hermitian');
    return
end

% Eigen-decomposition and compute square root
[V, D] = eig(X_opt);
eigenvalues = real(diag(D));
eigenvalues = max(eigenvalues, 0);  % Ensuring eigenvalues are non-negative
sqrt_X = V * diag(sqrt(eigenvalues));  % Square root matrix

% Initialize best solution tracking
best_power = inf;
W_best = [];

% Perform Gaussian randomizations
for t = 1:num_randomizations
    % Generate random complex vector z in \parial B_1(0)
    z = (randn(num_antennas, 1) + 1j * randn(num_antennas, 1));
    z = z / norm(z);

    % Generate candidate beamformer
    w_candidate = sqrt_X * z;
    
    % Compute received power at each user: p_k = |h_k^H * w|^2
    pk = zeros(num_users, 1);
    for k = 1:num_users
        pk(k) = abs(H(:, k)' * w_candidate)^2;
    end
    
    % Find scaling factor to satisfy all SNR constraints
    % We need: pk_scaled(k) >= rhs(k) for all k
    % pk_scaled(k) = alpha^2 * pk(k)
    % So: alpha^2 >= rhs(k) / pk(k) for all k
    alpha2 = max(rhs ./ max(pk, eps));  % Avoid division by zero
    
    % Scale the candidate
    w_scaled = sqrt(alpha2) * w_candidate;
    
    % Compute transmit power
    transmit_power = norm(w_scaled)^2;
    
    % Keep the solution with minimum power
    if transmit_power < best_power && isfinite(transmit_power)
        best_power = transmit_power;
        W_best = w_scaled;
    end
end

% Check if we found a valid solution
if isempty(W_best) || ~isfinite(best_power)
    % Fallback: use dominant eigenvector with scaling
    [~, max_idx] = max(eigenvalues);
    w_candidate = sqrt(eigenvalues(max_idx)) * V(:, max_idx);
    pk = zeros(num_users, 1);
    for k = 1:num_users
        pk(k) = abs(H(:, k)' * w_candidate)^2;
    end
    alpha2 = max(rhs ./ max(pk, eps));
    W_best = sqrt(alpha2) * w_candidate;
    best_power = norm(W_best)^2;
    metrics.status_message = 'Used dominant eigenvector (randomization failed)';
else
    metrics.status_message = 'Solved via SDR + Gaussian randomization';
end


W = W_best;

%% Compute final metrics using the chosen W
% Use centralized metrics computation for consistency
metrics = compute_beamformer_metrics(W, H, config);

% Add SDR-specific metrics
metrics.cvx_optval = cvx_optval;
metrics.cvx_optval_db = 10*log10(cvx_optval);
metrics.num_randomizations = num_randomizations;
metrics.converged = true;


%{
  Need SDR bound normalization if P_tr is provided in config for cvx optval
  because cvx_optval corresponds to the minimum power at which the QoS
%}

% Normalize SNR if P_tr is provided
if isfield(config, 'P_tr') && ~isempty(config.P_tr) && config.P_tr > 0
    % Scale the SDR bound: cvx_optval represents minimum power
    % If we normalize to P_tr, the bound scales proportionally
    X_opt_normalized = X_opt * config.P_tr / cvx_optval;
    SNR_cvx = zeros(num_users, 1);
    for k = 1:num_users
        % Compute SNR with normalized power
        SNR_cvx(k) = trace(Q{k} * X_opt_normalized) / noise_power(k);
    end
    SNR_cvx = real(SNR_cvx);
    min_SNR_cvx = min(SNR_cvx);
    metrics.cvx_min_snr = min_SNR_cvx;
    metrics.cvx_min_snr_db = 10*log10(metrics.cvx_min_snr);
end

% Add algorithm-specific status message
if ~metrics.feasible
    metrics.status_message = 'WARNING: QoS constraints not satisfied at actual power (SDR + randomization)';
end

end
