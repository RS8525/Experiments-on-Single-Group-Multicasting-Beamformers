function satisfied = check_power_constraint(W, P_max, tolerance)
% CHECK_POWER_CONSTRAINT Verify power constraint is satisfied
%
% Inputs:
%   W         - [M x N] beamforming matrix
%   P_max     - Maximum transmit power
%   tolerance - (optional) Tolerance for constraint (default = 1e-6)
%
% Outputs:
%   satisfied - true if ||W||_F^2 <= P_max + tolerance, false otherwise
%
% Example:
%   satisfied = check_power_constraint(W, 1.0);

if nargin < 3
    tolerance = 1e-6;
end

power = norm(W, 'fro')^2;
satisfied = (power <= P_max + tolerance);

end
