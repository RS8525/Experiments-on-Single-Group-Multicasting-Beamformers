function alpha = compute_scaling_factor(w, H, gamma, sigma_k_squared)
% COMPUTE_SCALING_FACTOR Computes the scaling factor to satisfy all SNR constraints
%
% Finds the maximum scaling factor alpha^k needed across all users k such
% that all SNR constraints are satisfied. If the maximum is less than 1,
% no scaling is needed.
%
% Inputs:
%   w              : Candidate precoding vector
%   H              : [num_antennas x num_users] channel matrix
%   gamma          : [num_users x 1] QoS/SNR targets
%   sigma_k_squared: [num_users x 1] noise power
%
% Output:
%   alpha : Scaling factor (>= 1 if scaling needed, otherwise 1)
%
% Formula: alpha = max_k { gamma_k * sigma_k^2 / |h_k^H * w|^2 }

% Compute received power at each user
recv_power = abs(H' * w).^2;  % [num_users x 1]

% Avoid division by zero
recv_power = max(recv_power, 1e-12);

% Compute required scaling factors for each user
alpha_k = (gamma .* sigma_k_squared) ./ recv_power;  % [num_users x 1]

% Take the maximum (worst-case user)
alpha = max(alpha_k);

%{
 % If alpha < 1, no scaling is needed
alpha = max(alpha, 1);
 
%}
end
