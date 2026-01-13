function set_global_seed(seed)
% SET_GLOBAL_SEED Set global random seed for reproducibility
%
% Inputs:
%   seed - Integer seed value for random number generators
%
% This function sets the seed for both MATLAB's default random number
% generator (rng) to ensure reproducible results across experiments.
%
% Example:
%   set_global_seed(42);
%   H = randn(4, 2);  % Will produce same result with same seed

validateattributes(seed, {'numeric'}, {'scalar', 'integer', 'nonnegative'});

% Set MATLAB's random number generator
rng(seed, 'twister');

fprintf('Global random seed set to: %d\n', seed);

end
