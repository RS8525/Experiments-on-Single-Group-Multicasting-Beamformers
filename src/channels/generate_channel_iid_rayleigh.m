function H = generate_channel_iid_rayleigh(num_antennas, num_users, num_streams, varargin)
% GENERATE_CHANNEL_IID_RAYLEIGH Generate IID Rayleigh fading channel
%
% Inputs:
%   num_antennas - Number of transmit antennas
%   num_users    - Number of users
%   num_streams  - Number of streams (optional, default = 1)
%   varargin - Optional name-value pairs:
%              'Normalize', true/false - Normalize channel (default: false)
%
% Outputs:
%   H - [num_antennas x num_users x num_streams] channel matrix
%       H(:, k, n) is the channel vector from transmitter to user k for stream n
%
% Example:
%   H = generate_channel_iid_rayleigh(8, 4);        % [8 x 4 x 1]
%   H = generate_channel_iid_rayleigh(8, 4, 2);     % [8 x 4 x 2]
%   H = generate_channel_iid_rayleigh(8, 4, 1, 'Normalize', true);

% Parse inputs
if nargin < 3
    num_streams = 1;
end

p = inputParser;
addParameter(p, 'Normalize', false, @islogical);
parse(p, varargin{:});
normalize = p.Results.Normalize;

% Generate complex Gaussian channel (Rayleigh fading)
% Each element ~ CN(0, 1)
H = (randn(num_antennas, num_users, num_streams) + 1j * randn(num_antennas, num_users, num_streams)) / sqrt(2);

% Optional normalization
if normalize
    for n = 1:num_streams
        for k = 1:num_users
            H(:, k, n) = H(:, k, n) / norm(H(:, k, n));
        end
    end
end

end
