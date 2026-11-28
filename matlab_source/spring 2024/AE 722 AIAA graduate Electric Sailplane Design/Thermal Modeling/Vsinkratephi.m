function V_scphi = Vsinkratephi(C_Dmax, C_Lmax, r, g, phi, V_k, rho, W2S)

    for i = 1:length(phi)
        V_scphi(i,:) = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S) .* 1./(r .* g.* (sin(phi(1,i))./(V_k.^2))); 
    end
end