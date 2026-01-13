function snr = compute_snr(W, H, config)
% COMPUTE_SNR Compute SNR for each user
%
% Inputs:
%   W      - [M x N] beamforming matrix (M antennas, N streams)
%   H      - [M x K] or [M x K x N] channel matrix (K users)
%   config - struct with fields:
%            .noise_power : noise power (default = 1.0)
%
% Outputs:
%   snr - [K x 1] vector of SNR values for each user
%
% Example:
%   snr = compute_snr(W, H, config);
%   snr_db = 10*log10(snr);

% Extract dimensions
[M, N] = size(W);
if ndims(H) == 2
    [~, K] = size(H);
    H = reshape(H, M, K, 1);
else
    [~, K, ~] = size(H);
end

% Extract noise power
if isfield(config, 'noise_power')
    sigma2 = config.noise_power;
else
    sigma2 = 1.0;
end

% Compute SNR for each user
snr = zeros(K, 1);
for k = 1:K
    signal_power = 0;
    for n = 1:N
        signal_power = signal_power + abs(H(:, k, n)' * W(:, n))^2;
    end
    snr(k) = signal_power / sigma2;
end

end
