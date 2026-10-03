function V_sc_phi = Vsinkrate_phi(C_Dmax, C_Lmax, phi, rho, W2S)
V_sc_phi_W2S1 = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,1)); % Sink rate m/s at specific roll angles at wing loading 1
V_sc_phi_W2S2 = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,2));
V_sc_phi_W2S3 = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,3));
V_sc_phi_W2S4 = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,4));
V_sc_phi_W2S5 = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,5));
V_sc_phi_W2S6 = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,6));
V_sc_phi_W2S7 = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,7));
V_sc_phi = [V_sc_phi_W2S1; V_sc_phi_W2S2; V_sc_phi_W2S3; V_sc_phi_W2S4; V_sc_phi_W2S5; V_sc_phi_W2S6; V_sc_phi_W2S7];

end