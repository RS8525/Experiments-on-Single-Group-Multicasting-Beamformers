function H = generate_channel_iid_rayleigh(M, K, N, varargin)
% GENERATE_CHANNEL_IID_RAYLEIGH Generate IID Rayleigh fading channel
%
% Inputs:
%   M - Number of transmit antennas
%   K - Number of users
%   N - Number of streams (optional, default = 1)
%   varargin - Optional name-value pairs:
%              'Normalize', true/false - Normalize channel (default: false)
%
% Outputs:
%   H - [M x K x N] channel matrix
%       H(:, k, n) is the channel vector from transmitter to user k for stream n
%
% Example:
%   H = generate_channel_iid_rayleigh(8, 4);        % [8 x 4 x 1]
%   H = generate_channel_iid_rayleigh(8, 4, 2);     % [8 x 4 x 2]
%   H = generate_channel_iid_rayleigh(8, 4, 1, 'Normalize', true);

% Parse inputs
if nargin < 3
    N = 1;
end

p = inputParser;
addParameter(p, 'Normalize', false, @islogical);
parse(p, varargin{:});
normalize = p.Results.Normalize;

% Generate complex Gaussian channel (Rayleigh fading)
% Each element ~ CN(0, 1)
H = (randn(M, K, N) + 1j * randn(M, K, N)) / sqrt(2);

% Optional normalization
if normalize
    for n = 1:N
        for k = 1:K
            H(:, k, n) = H(:, k, n) / norm(H(:, k, n));
        end
    end
end

end
