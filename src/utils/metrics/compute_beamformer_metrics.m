function metrics = compute_beamformer_metrics(W, H, config)
% COMPUTE_BEAMFORMER_METRICS Compute standardized performance metrics for beamforming
%
% This function provides centralized metric computation for all beamforming
% algorithms. It optionally normalizes SNR metrics to a reference power budget
% P_tr for fair comparison across algorithms, while keeping the actual
% beamforming vector W and its power unchanged.
%
% Inputs:
%   W      - [num_antennas x 1] beamforming vector
%   H      - [num_antennas x num_users] channel matrix
%   config - struct with fields:
%            .noise_power : [num_users x 1] noise power (required)
%            .gamma       : [num_users x 1] QoS/SNR targets (required)
%            .P_tr        : Power budget for SNR normalization (optional)
%                           If provided, SNR metrics computed as if W scaled to P_tr
%                           If omitted, SNR computed from actual W
%
% Outputs:
%   metrics - struct with fields:
%             .final_power     : actual transmit power ||W||^2
%             .snr             : [num_users x 1] SNR per user
%             .min_snr         : minimum SNR across users
%             .snr_db          : [num_users x 1] SNR in dB
%             .min_snr_db      : minimum SNR in dB
%             .power_db        : transmit power in dB
%             .rate            : sum rate in bits/s/Hz
%             .feasible        : true if actual SNR satisfies QoS constraints
%
% Note on Normalization:
%   - When config.P_tr is provided, SNR metrics reflect performance at P_tr
%   - Feasibility is always checked against ACTUAL SNR (from actual W)
%   - final_power always reflects actual power of W
%   - This allows fair SNR comparison across algorithms at same power budget
%
% Example:
%   % Without normalization (standard metrics)
%   config = struct('noise_power', ones(4,1), 'gamma', ones(4,1));
%   metrics = compute_beamformer_metrics(W, H, config);
%
%   % With normalization (SNR at 1W for fair comparison)
%   config.P_tr = 1.0;
%   metrics = compute_beamformer_metrics(W, H, config);

%% Input validation
[num_antennas, num_users] = size(H);

assert(isvector(W) && length(W) == num_antennas, ...
    'W must be [num_antennas x 1] vector');

assert(isfield(config, 'noise_power'), ...
    'config.noise_power is required');
sigma_k_squared = config.noise_power(:);
assert(length(sigma_k_squared) == num_users, ...
    'noise_power must have num_users elements');

assert(isfield(config, 'gamma'), ...
    'config.gamma is required');
gamma = config.gamma(:);
assert(length(gamma) == num_users, ...
    'gamma must have num_users elements');

% Validate P_tr if provided
if isfield(config, 'P_tr') && ~isempty(config.P_tr)
    assert(isscalar(config.P_tr) && config.P_tr > 0 && isfinite(config.P_tr), ...
        'P_tr must be a positive finite scalar');
end

%% Compute actual power
W = W(:);  % Ensure column vector
P_actual = norm(W)^2;
metrics.final_power = P_actual;
metrics.power_db = 10*log10(P_actual);

%% Determine beamformer for SNR computation
if isfield(config, 'P_tr') && ~isempty(config.P_tr) && config.P_tr > 0
    % Normalize to reference power for fair SNR comparison
    % W_normalized = W * sqrt(P_tr / P_actual)
    scale_factor = sqrt(config.P_tr / P_actual);
    W_metric = W * scale_factor;
else
    % Use actual beamformer
    W_metric = W;
end

%% Compute SNR metrics (from W_metric)
recv_power_metric = abs(H' * W_metric).^2;  % [num_users x 1]
snr_metric = recv_power_metric ./ sigma_k_squared;

metrics.snr = snr_metric;
metrics.min_snr = min(snr_metric);
metrics.snr_db = 10*log10(snr_metric);
metrics.min_snr_db = 10*log10(metrics.min_snr);

%% Compute rate (from normalized SNR)
metrics.rate = sum(log2(1 + snr_metric));

%% Check feasibility against ACTUAL SNR (not normalized)
% This flags beamformers that violate QoS at their actual power
recv_power_actual = abs(H' * W).^2;
snr_actual = recv_power_actual ./ sigma_k_squared;

tol_feas = 1e-6;
metrics.feasible = all(snr_actual >= gamma * (1 - tol_feas));

end
