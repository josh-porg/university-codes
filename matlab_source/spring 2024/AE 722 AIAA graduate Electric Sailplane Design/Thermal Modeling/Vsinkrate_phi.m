function V_sc_phi = Vsinkrate_phi(C_Dmax, C_Lmax, phi, rho, W2S)

for i = 1:length(W2S)
V_sc_phi(i,:) = (C_Dmax ./ (C_Lmax^(3/2) .* (cos(phi)).^(3/2))) .* sqrt((2 ./ rho) .* W2S(1,i)); % Sink rate m/s at specific roll angles at wing loading 1
end

end