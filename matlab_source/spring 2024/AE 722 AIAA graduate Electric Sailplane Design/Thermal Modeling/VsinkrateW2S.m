function V_scW2S = VsinkrateW2S(C_Dmax, C_Lmax, r, g, phi, V_k, rho, W2S)

    for i = 1:length(W2S)
        V_scW2S(i,:) = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,i)) .* 1./(r .* g.* (sin(phi)./(V_k.^2))); 
    end
end