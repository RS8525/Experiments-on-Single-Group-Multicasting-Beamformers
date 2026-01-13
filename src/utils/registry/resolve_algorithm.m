function algo_fn = resolve_algorithm(name)
% RESOLVE_ALGORITHM Map algorithm name to function handle.
%
% Inputs:
%   name - Algorithm name (string)
%
% Outputs:
%   algo_fn - Function handle

if isempty(name)
    error('Algorithm name is required');
end

name = char(lower(string(name)));

switch name
    case 'sdr_beamformer'
        fn_name = 'sdr_beamformer';
    case 'ff_c2'
        fn_name = 'ff_c2';
    case 'rc_c2'
        fn_name = 'rc_c2';
    case 'sbfc'
        fn_name = 'sbfc';
    case 'rc_c2_update'
        fn_name = 'rc_c2_update';
    otherwise
        error('Unknown algorithm: %s', name);
end

if exist(fn_name, 'file') ~= 2
    error('Algorithm not implemented: %s', fn_name);
end

algo_fn = str2func(fn_name);

end
