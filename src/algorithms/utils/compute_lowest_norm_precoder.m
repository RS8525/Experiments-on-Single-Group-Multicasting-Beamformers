function w = compute_lowest_norm_precoder(h_i, h_j, gamma_i, gamma_j, sigma_i_sq, sigma_j_sq)
% COMPUTE_LOWEST_NORM_PRECODER Computes the minimum norm precoder for two users
%
% This function computes the precoding vector w^{i,j} that satisfies the
% SNR constraints for users i and j with the lowest transmit power.
% The specific formula depends on the relationship between the channels.
%
% Inputs:
%   h_i, h_j     : Channel vectors for users i and j
%   gamma_i, gamma_j : QoS/SNR targets for users i and j
%   sigma_i_sq, sigma_j_sq : Noise power for users i and j
%
% Output:
%   w : Precoding vector that meets SNR constraints with minimum norm
%
% Reference: Equations (12), (14), (17), (19) from the paper

% Precompute constants
c_i = gamma_i * sigma_i_sq;
c_j = gamma_j * sigma_j_sq;

% Compute channel properties
norm_h_i_sq = norm(h_i)^2;
norm_h_j_sq = norm(h_j)^2;
inner_product = h_j' * h_i;
inner_product_sq = abs(inner_product)^2;

% Case distinction based on which constraints are tight
% Case 1 (Eq. 12): Only constraint i is tight
% Condition (Eq. 15): |h_j^H h_i|^2 <= (c_j/c_i) * ||h_i||^4
if inner_product_sq >= (c_j / c_i) * norm_h_i_sq^2
    % Equation (12): w = sqrt(c_i) / ||h_i|| * h_i
    w = sqrt(c_i / norm_h_i_sq) * h_i;
    
% Case 2 (Eq. 14): Only constraint j is tight  
% Condition (Eq. 13): |h_i^H h_j|^2 <= (c_i/c_j) * ||h_j||^4
elseif inner_product_sq >= (c_i / c_j) * norm_h_j_sq^2
    % Equation (14): w = sqrt(c_j) / ||h_j|| * h_j
    w = sqrt(c_j / norm_h_j_sq) * h_j;
    
% Case 3 (Eq. 17): Both constraints are tight (collinear case)
else
    % Form H matrix: H = [h_i'; h_j'] (2 x M matrix)
    H = [h_i'; h_j'];
    
    % Compute Gram matrix: G = H * H^H (2 x 2 matrix)
    G = H * H';
    
    % Invert Gram matrix
    G_inv = inv(G);
    
    % Phase optimization (Eq. 19):
    % To minimize transmit power, set phi_i = 1 and optimize phi_j
    % Optimal phase: phi_j = exp(j*(pi - angle(z_3)))
    % where z_3 = G_inv(1,2) is the off-diagonal element
    z_3 = G_inv(1, 2);
    phi_i = 1;
    phi_j = exp(1j * (pi - angle(z_3)));
    
    % Equation (17): w = H^H * (H*H^H)^(-1) * [sqrt(c_i)*phi_i; sqrt(c_j)*phi_j]
    w = H' * G_inv * [sqrt(c_i) * phi_i; sqrt(c_j) * phi_j];
end

end
