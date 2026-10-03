function [Pd] = Fuze(V, r_tgt_missile, L, K, S, phi, time_elapsed, x0_t, y0_t, Vt_mag, T_pos)
    % This function calculates the probability of damage (Pd) based on the current 
    % interceptor state and target parameters using an exponential model.
    
    % Compute target's current position.
    St = T_pos;
    % Relative target position with respect to the interceptor.
    St_prime = St - S;
    Range = norm(St_prime);
    
    % Compute the off-axis angle from the interceptor velocity to target velocity.
    Vt = [-Vt_mag 0]; 
    top = cross([V 0], [Vt 0]);
    sin_theta_off = top(3)/(norm(V)*norm(Vt));
    cos_psi = sin_theta_off;
    A = 2 * r_tgt_missile * L; 
    Ae = A * cos_psi;
    if Ae < pi * r_tgt_missile^2
        Ae = pi * r_tgt_missile^2;
    end
    
    % Compute the effective spreading area at the target range.
    As = 2 * pi * (Range^2) * (1 - cos(phi));
    % Exponential fuze/damage probability equation.
    Pd = 1 - exp(-K * (Ae/As));
end