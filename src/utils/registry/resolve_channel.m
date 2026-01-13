function channel_fn = resolve_channel(channel_config)
% RESOLVE_CHANNEL Map channel config to a generator function handle.
%
% Inputs:
%   channel_config - Struct with channel settings
%
% Outputs:
%   channel_fn - Function handle (M, K, N) -> H

if nargin < 1 || isempty(channel_config)
    channel_config = struct();
end

channel_type = 'iid_rayleigh';
if isfield(channel_config, 'type')
    channel_type = channel_config.type;
end

normalize = false;
if isfield(channel_config, 'normalize')
    normalize = logical(channel_config.normalize);
end

channel_type = char(lower(string(channel_type)));

switch channel_type
    case 'iid_rayleigh'
        if exist('generate_channel_iid_rayleigh', 'file') ~= 2
            error('Channel generator not found: generate_channel_iid_rayleigh');
        end
        channel_fn = @(M, K, N) generate_channel_iid_rayleigh(M, K, N, 'Normalize', normalize);
    otherwise
        error('Unknown channel type: %s', channel_type);
end

end
