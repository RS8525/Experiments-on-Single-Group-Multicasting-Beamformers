function [w, h_orth] = sbfc(h, c, gamma, N, rc_c2)
    % SBFC Successive Beamforming Filter Computation
    %
    % Syntax:
    %   [w, h_orth] = sbfc(h, c, gamma, N, rc_c2)
    %
    % Input:
    %   H      - [num_antennas x num_users] channel matrix
    %   c       - Cost vector (K x 1) for user selection
    %   gamma   - SNR threshold for exit condition
    %   N       - Maximum iterations
    %   rc_c2   - Reference constraint computation function handle
    %
    % Output:
    %   w       - Beamforming filter weights (M x 1)
    %   h_orth  - Orthogonalized channels (M x K)
    
    % Step 1: Save original channels
    h_orig = H;
    [num_antennas, num_users] = size(H);
    h_orth = h_orig;
    
    % Step 2: Select weakest user (most violated constraint)
    [~, l1] = min(c .* vecnorm(h_orth, 2, 1).^2);
    
    % Step 3: Set initial filter
    w = (h_orth(:, l1) / norm(h_orth(:, l1), 2));
    
    % Step 4: Initialize active user set
    U = setdiff(1:num_users, l1);
    
    % Main iteration loop
    for n = 2:N
        % Step 6: Find most violated constraint
        snr_k = zeros(1, num_users);
        for k = U
            snr_k(k) = abs(w' * h_orth(:, k))^2;
        end
        
        snr_all = abs(w' * h_orth).^2;
        snr_ratio = snr_k(U) ./ snr_all(U);
        
        [~, idx] = max(snr_ratio);
        l_n = U(idx);
        
        % Step 7: Exit if all constraints met
        if snr_ratio(idx) > 1
            break;
        end
        
        % Step 8: Orthogonalization
        for k = U
            h_orth(:, k) = rc_c2(h_orth(:, k), h_orth(:, l_n));
        end
        
        % Step 9: Update filter
        alpha_n = compute_alpha(snr_k(l_n), snr_all(l_n), w, h_orth(:, l_n));
        w = w + alpha_n * (h_orth(:, l_n) / norm(h_orth(:, l_n), 2));
        
        % Step 10: Update active set
        U = setdiff(U, l_n);
    end
end

function alpha = compute_alpha(snr_k, snr_l, w, h_l)
    % Compute step size for filter update (equations 27-28)
    h_l_norm = norm(h_l, 2);
    alpha = (sqrt(snr_k) - sqrt(snr_l)) / h_l_norm;
end