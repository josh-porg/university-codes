function Pd = estimatePd(Active, n, phi_deg, delta_x_d, f, r_tgt_missile, L_target, K)
    % This function estimates Pd based on detector information and given phi
    % Zaeem Siddiqui
    
    % Convert phi from degrees to radians
    phi = deg2rad(phi_deg);   
    
    % Estimate image width on detector array (in meters)
    x_illum = (Active - 1) * delta_x_d / 2;
    
    % Estimate range
    Le = L_target;  % Assume full target length (most conservative case, perpendicular)
    r_estimated = (Le * f) / x_illum;
    
    % Estimate effective area
    Ae = 2 * r_tgt_missile * L_target;    % Side-on projected area
    if Ae < pi * r_tgt_missile^2
        Ae = pi * r_tgt_missile^2;  % Use minimum value of the projected area
    end
    
    % Fragment cone surface area at estimated range
    Omega = 2 * pi * (1 - cos(phi));
    As = Omega * r_estimated^2;
    
    % Calculate Pd estimate
    Pd = 1 - exp(-K * (Ae / As));
    
    % Display results (optional)
    fprintf('--- Estimated Pd ---\n');
    fprintf('Center Detector Index (n): %d\n', n);
    fprintf('Detectors Lit (Active): %d\n', Active);
    fprintf('Estimated Range: %.3f m\n', r_estimated);
    fprintf('Estimated Pd: %.3f\n', Pd);
end
