function rate = compute_sum_rate(snr)
% COMPUTE_SUM_RATE Compute sum rate from SNR values
%
% Inputs:
%   snr - [K x 1] vector of SNR values (linear, not dB)
%
% Outputs:
%   rate - Sum rate in bits/s/Hz
%
% Formula:
%   rate = sum_k log2(1 + snr_k)
%
% Example:
%   snr = [10; 15; 20];  % Linear SNR values
%   rate = compute_sum_rate(snr);

rate = sum(log2(1 + snr));

end
