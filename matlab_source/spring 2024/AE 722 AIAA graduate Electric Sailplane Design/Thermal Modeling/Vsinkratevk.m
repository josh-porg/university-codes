function V_scvk = Vsinkratevk(C_Dmax, C_Lmax, r, g, phi, V_k, rho, W2S)

    for i = 1:length(V_k)
        V_scvk(i,:) = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S) .* 1./(r .* g.* (sin(phi)./(V_k(1,i).^2)));
    end 
end