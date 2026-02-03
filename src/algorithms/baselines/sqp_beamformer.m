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


%% Initial guess

% Use RC-C2 as initialization
w0 = rc_c2(H, config);


% Convert to real representation for fmincon
% w = w_real + 1j*w_imag -> x = [w_real; w_imag]
x0 = [real(w0); imag(w0)];  % [2*num_antennas x 1]

%% Define optimization problem

% Objective function: minimize ||w||^2 with gradient
objective = @(x) sqp_objective(x);

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



% Determine convergence status based on exitflag
% exitflag > 0: converged successfully
% exitflag = 0: max iterations reached
% exitflag < 0: failed


% Status message based on exitflag
if exitflag > 0
    if ~metrics.feasible
        metrics.status_message = 'SQP converged, WARNING: QoS constraints not satisfied';
    else
        metrics.status_message = sprintf('SQP converged (%d iterations)', output.iterations);
    end
elseif exitflag == 0
    metrics.status_message = sprintf('SQP: max iterations reached (%d)', max_iter);
elseif exitflag < 0
    metrics = compute_beamformer_metrics(w0, H, config);
    W = w0;
    metrics.status_message = 'SQP: failed falling back to rc_c2 initialization';
else
    metrics.status_message = sprintf('SQP failed (exitflag=%d)', exitflag);
end

% Add SQP-specific fields
metrics.converged = (exitflag > 0);
metrics.solve_time = solve_time;
metrics.iterations = output.iterations;
metrics.exitflag = exitflag;
metrics.fval = fval;

end

%% Helper function: Objective with gradient
function [f, g] = sqp_objective(x)
% SQP_OBJECTIVE Objective function for SQP: minimize ||w||^2
%
% Inputs:
%   x - [2*M x 1] real representation [w_real; w_imag]
%
% Outputs:
%   f - objective value ||w||^2
%   g - gradient [2*M x 1]

% Objective: f = ||w||^2 = w_real'*w_real + w_imag'*w_imag
f = norm(x)^2;

% Gradient: df/dx = 2*x
g = 2 * x;

end

%% Helper function: Nonlinear constraints with gradients
function [c, ceq, gc, gceq] = sqp_constraints(x, H, gamma, sigma2)
% SQP_CONSTRAINTS Nonlinear constraints for SQP beamformer
%
% Inputs:
%   x              - [2*M x 1] real representation [w_real; w_imag]
%   H              - [M x K] channel matrix
%   gamma          - [K x 1] QoS targets
%   sigma2         - [K x 1] noise powers
%
% Outputs:
%   c   - [K x 1] inequality constraints (c <= 0)
%   ceq - [] no equality constraints
%   gc  - [2*M x K] gradient of inequality constraints
%   gceq- [] no equality constraint gradients

[num_antennas, num_users] = size(H);

w = x(1:num_antennas) + 1j*x(num_antennas+1:end);

s = H' * w;                 % Kx1, s_k = h_k^H w
recv = abs(s).^2;           % Kx1

c   = gamma(:).*sigma2(:) - recv;
ceq = [];

% Gradient of recv_k = |s_k|^2 wrt [Re(w); Im(w)] is:
% 2*[Re(h_k*s_k); Im(h_k*s_k)]
G = 2 * H .* (s.' );        % MxK, column k is 2*h_k*s_k   (NOTE: no conj)

gc = -[real(G); imag(G)];   % (2M)xK because c = ... - recv
gceq = [];

% Defensive: ensure real outputs
c  = real(c);
gc = real(gc);

end
