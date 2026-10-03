function [C_Lopt,V_c,C_D0] = CLoptimum(L_Dmax, A, e, rho, r, g, W2S)

C_Lmax = (pi .* A .* e) ./ (2 .* L_Dmax); % Maximum lift coefficient
C_Dmax = 2 .* (C_Lmax.^2)./(pi .* A .* e); % Maximum drag coefficient
C_D0 = C_Dmax ./ 2; % Parasite drag coefficient

V_T = HorstmannThermal(r);
V_sc = Vsinkrate(C_Dmax, C_Lmax, r, g, rho, W2S);
[V_ca1, V_ca2, V_cb1, V_cb2] = Vclimbspeed(V_T, V_sc);
V_c = [V_ca1; V_ca2; V_cb1; V_cb2];

syms C_Lopt 

eqn = 0 == C_D0(1,1) - (C_Lopt^2) ./ (pi.*A(1,1).*e) - (V_c(:,:) .* C_Lopt^(3/2)) ./ (2 .* sqrt((2/rho) .* W2S(1,1)));
for j = 1:4  
    for i = 1:length(V_ca1) 
         C_L(j,i) = double(vpasolve(eqn(j,i),C_Lopt));
    end
end

end