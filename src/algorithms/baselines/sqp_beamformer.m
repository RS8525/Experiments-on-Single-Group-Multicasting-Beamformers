function [W, metrics] = sqp_beamformer(H, config)
% SQP_BEAMFORMER Sequential Quadratic Programming baseline beamformer
%
% Solves the multicast beamforming problem using MATLAB's fmincon with
% Sequential Quadratic Programming (SQP) algorithm. This is a standard
% nonlinear optimization baseline.
%
% Inputs:
%   H      - [num_antennas x num_users] channel matrix
%            H(:,k) is the channel vector from transmitter to user k
%   config - struct with fields:
%            .gamma           : [num_users x 1] QoS/SNR targets (required)
%            .noise_power     : [num_users x 1] noise power (required)
%            .sqp_max_iter    : maximum iterations (optional, default = 1000)
%            .sqp_tol         : tolerance (optional, default = 1e-6)
%            .sqp_display     : 'off', 'iter', 'final' (optional, default = 'off')
%            .seed            : random seed for initialization (optional)
%
% Outputs:
%   W       - [num_antennas x 1] beamforming vector
%   metrics - struct with fields:
%             .final_power     : final transmit power
%             .snr             : [num_users x 1] SNR per user
%             .min_snr         : minimum SNR across users
%             .rate            : sum rate in bits/s/Hz
%             .feasible        : true if QoS constraints satisfied
%             .converged       : true if optimizer converged
%             .solve_time      : solver time in seconds
%             .iterations      : number of iterations
%             .exitflag        : fmincon exit flag
%             .status_message  : descriptive status string
%
% Requirements:
%   MATLAB Optimization Toolbox
%
% Method:
%   Solve: minimize ||w||^2
%          subject to: |h_k' * w|^2 >= gamma_k * sigma_k^2, for all k
%
%   Uses fmincon with 'sqp' algorithm to solve this constrained optimization.
%
% Reference:
%   SQP is a standard nonlinear optimization method for constrained problems

%% Configuration and parameter extraction
[num_antennas, num_users] = size(H);

% Extract QoS requirements
if ~isfield(config, 'gamma')
    error('sqp_beamformer:missingParameter', 'config.gamma is required');
end
gamma = config.gamma(:);
assert(length(gamma) == num_users, ...
    'gamma must be [num_users x 1] vector, got [%d x 1]', length(gamma));

% Extract noise power
if ~isfield(config, 'noise_power')
    error('sqp_beamformer:missingParameter', 'config.noise_power is required');
end
sigma_k_squared = config.noise_power(:);
assert(length(sigma_k_squared) == num_users, ...
    'noise_power must be [num_users x 1] vector, got [%d x 1]', length(sigma_k_squared));

% SQP parameters with defaults
if isfield(config, 'sqp_max_iter')
    max_iter = config.sqp_max_iter;
else
    max_iter = 1000;
end

if isfield(config, 'sqp_tol')
    tol = config.sqp_tol;
else
    tol = 1e-6;
end

if isfield(config, 'sqp_display')
    display_mode = config.sqp_display;
else
    display_mode = 'off';
end

% Set random seed for initialization if provided
if isfield(config, 'seed')
    if exist('set_global_seed', 'file')
        set_global_seed(config.seed);
    else
        rng(config.seed, 'twister');
    end
end

%% Initial guess

% Use MRT (Maximum Ratio Transmission) direction as initialization
% MRT: w = sum(h_k) / ||sum(h_k)||
h_sum = sum(H, 2);  % Sum across all user channels
w0 = h_sum / norm(h_sum);

% Scale to roughly satisfy average QoS (rough initialization)
avg_gamma = mean(gamma);
avg_noise = mean(sigma_k_squared);
recv_power_init = abs(H' * w0).^2;
avg_recv_power = mean(recv_power_init);
scale_init = sqrt(avg_gamma * avg_noise / max(avg_recv_power, eps));
w0 = w0 * scale_init;

% Convert to real representation for fmincon
% w = w_real + 1j*w_imag -> x = [w_real; w_imag]
x0 = [real(w0); imag(w0)];  % [2*num_antennas x 1]

%% Define optimization problem

% Objective function: minimize ||w||^2
objective = @(x) norm(x)^2;

% Nonlinear constraints: |h_k' * w|^2 >= gamma_k * sigma_k^2
% In fmincon format: c(x) <= 0, so we use:
% gamma_k * sigma_k^2 - |h_k' * w|^2 <= 0  (flip inequality)
nonlcon = @(x) sqp_constraints(x, H, gamma, sigma_k_squared);

% No linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% No bounds
lb = [];
ub = [];

%% Configure fmincon options
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', max_iter, ...
    'MaxFunctionEvaluations', max_iter * 100, ...
    'OptimalityTolerance', tol, ...
    'ConstraintTolerance', tol, ...
    'StepTolerance', tol, ...
    'Display', display_mode, ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true);

%% Solve with fmincon
tic;
[x_opt, fval, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
solve_time = toc;

%% Extract solution
% Convert back to complex beamformer
w_real = x_opt(1:num_antennas);
w_imag = x_opt(num_antennas+1:end);
W = w_real + 1j * w_imag;

%% Compute metrics
metrics = compute_beamformer_metrics(W, H, config);

% Add SQP-specific fields
metrics.solve_time = solve_time;
metrics.iterations = output.iterations;
metrics.exitflag = exitflag;
metrics.fval = fval;

% Determine convergence status based on exitflag
% exitflag > 0: converged successfully
% exitflag = 0: max iterations reached
% exitflag < 0: failed
metrics.converged = (exitflag > 0);

% Status message based on exitflag
if exitflag > 0
    if ~metrics.feasible
        metrics.status_message = 'SQP converged, WARNING: QoS constraints not satisfied';
    else
        metrics.status_message = sprintf('SQP converged (%d iterations)', output.iterations);
    end
elseif exitflag == 0
    metrics.status_message = sprintf('SQP: max iterations reached (%d)', max_iter);
elseif exitflag == -2
    metrics.status_message = 'SQP: no feasible point found';
else
    metrics.status_message = sprintf('SQP failed (exitflag=%d)', exitflag);
end

end

%% Helper function: Nonlinear constraints with gradients
function [c, ceq, gc, gceq] = sqp_constraints(x, H, gamma, sigma_k_squared)
% SQP_CONSTRAINTS Nonlinear constraints for SQP beamformer
%
% Inputs:
%   x              - [2*M x 1] real representation [w_real; w_imag]
%   H              - [M x K] channel matrix
%   gamma          - [K x 1] QoS targets
%   sigma_k_squared- [K x 1] noise powers
%
% Outputs:
%   c   - [K x 1] inequality constraints (c <= 0)
%   ceq - [] no equality constraints
%   gc  - [2*M x K] gradient of inequality constraints
%   gceq- [] no equality constraint gradients

[num_antennas, num_users] = size(H);

% Extract real and imaginary parts
w_real = x(1:num_antennas);
w_imag = x(num_antennas+1:end);
w = w_real + 1j * w_imag;

% Compute constraints: c_k = gamma_k * sigma_k^2 - |h_k' * w|^2 <= 0
% We want |h_k' * w|^2 >= gamma_k * sigma_k^2
c = zeros(num_users, 1);
gc = zeros(2*num_antennas, num_users);

for k = 1:num_users
    h_k = H(:, k);
    
    % Constraint value
    recv_power = abs(h_k' * w)^2;
    c(k) = gamma(k) * sigma_k_squared(k) - recv_power;
    
    % Gradient computation
    % d/dw_real |h' * w|^2 = 2 * real(h .* conj(h' * w))
    % d/dw_imag |h' * w|^2 = 2 * imag(h .* conj(h' * w))
    
    h_w = h_k' * w;
    grad_real = 2 * real(h_k * conj(h_w));
    grad_imag = 2 * imag(h_k * conj(h_w));
    
    % Since c_k = ... - |h_k' * w|^2, gradient has negative sign
    gc(1:num_antennas, k) = -grad_real;
    gc(num_antennas+1:end, k) = -grad_imag;
end

% No equality constraints
ceq = [];
gceq = [];

end
