function V_sc = Vsinkrate(C_Dmax, C_Lmax, r, g, rho, W2S)

for i = 1:length(W2S)
    V_sc(i,:) = (C_Dmax ./ C_Lmax^(3/2)) .* sqrt((2 ./ rho) .* W2S(1,i)).* 1./((sqrt(1 - (((2 ./ rho) .* W2S(1,i)) .* (1 ./ (r .* g .* C_Lmax))) .^2)).^(3/2)); % Sink rate m/s at various turn radii and wing loadings
end

end