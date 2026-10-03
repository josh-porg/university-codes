function V_sc = Vsinkrate(C_D_0, C_L, A, e, r, g, rho, W2S)
% 
% V_sc_W2S1 = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,1)) .* (1 - (((2 ./ rho) .* W2S(1,1)) .* (1 ./ (r .* g .* C_Lmax))).^2).^(-3/4); % Sink rate m/s at various turn radii at wing loading 1
% V_sc_W2S2 = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,2)) .* (1 - (((2 ./ rho) .* W2S(1,2)) .* (1 ./ (r .* g .* C_Lmax))).^2).^(-3/4); % Sink rate m/s at various turn radii at wing loading 2
% V_sc_W2S3 = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,3)) .* (1 - (((2 ./ rho) .* W2S(1,3)) .* (1 ./ (r .* g .* C_Lmax))).^2).^(-3/4); % Sink rate m/s at various turn radii at wing loading 3
% V_sc_W2S4 = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,4)) .* (1 - (((2 ./ rho) .* W2S(1,4)) .* (1 ./ (r .* g .* C_Lmax))).^2).^(-3/4); % Sink rate m/s at various turn radii at wing loading 4
% V_sc_W2S5 = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,5)) .* (1 - (((2 ./ rho) .* W2S(1,5)) .* (1 ./ (r .* g .* C_Lmax))).^2).^(-3/4); % Sink rate m/s at various turn radii at wing loading 5
% V_sc_W2S6 = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,6)) .* (1 - (((2 ./ rho) .* W2S(1,6)) .* (1 ./ (r .* g .* C_Lmax))).^2).^(-3/4); % Sink rate m/s at various turn radii at wing loading 6
% V_sc_W2S7 = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,7)) .* (1 - (((2 ./ rho) .* W2S(1,7)) .* (1 ./ (r .* g .* C_Lmax))).^2).^(-3/4); % Sink rate m/s at various turn radii at wing loading 7

V_sc_W2S1 = sinkRate(C_L, C_D_0, A, e, W2S(1,1), rho, r, g);
V_sc_W2S2 = sinkRate(C_L, C_D_0, A, e, W2S(1,2), rho, r, g);
V_sc_W2S3 = sinkRate(C_L, C_D_0, A, e, W2S(1,3), rho, r, g);
V_sc_W2S4 = sinkRate(C_L, C_D_0, A, e, W2S(1,4), rho, r, g);
V_sc_W2S5 = sinkRate(C_L, C_D_0, A, e, W2S(1,5), rho, r, g);
V_sc_W2S6 = sinkRate(C_L, C_D_0, A, e, W2S(1,6), rho, r, g);
V_sc_W2S7 = sinkRate(C_L, C_D_0, A, e, W2S(1,7), rho, r, g);


    function V_sink = sinkRate(C_L, C_D_0, A, e, W2S, rho, r, g)
        V_sink = (C_L.^2/(pi*A*e) + C_D_0) * C_L.^(-3/2) * sqrt(2*W2S/rho) * (1 - (2*W2S./(rho*r*g*C_L)).^2).^(3/4);
    end

V_sc = [V_sc_W2S1; V_sc_W2S2; V_sc_W2S3; V_sc_W2S4; V_sc_W2S5; V_sc_W2S6; V_sc_W2S7];

end