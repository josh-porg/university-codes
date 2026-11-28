function V_avg = Vavgfsd(V_c, C_D0, A, e, rho, W2S, C_L)

for i = 1:4
    V_avg= (V_c .* C_L(i,:).^(1/2)) ./ (C_D0 .* C_L(i,:).^(3/2) + (1./(pi.*A.*e)) .* C_L(i,:).^(1/2) + (V_c ./ sqrt((2/rho).*W2S)));
end
end
