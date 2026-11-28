function [V_ca1, V_ca2, V_cb1, V_cb2] = Vclimbspeed(V_T, V_sc)

 V_ca1 = V_T(1,:) - -V_sc(:,:);
 V_ca2 = V_T(2,:) - -V_sc(:,:);
 V_cb1 = V_T(3,:) - -V_sc(:,:);
 V_cb2 = V_T(4,:) - -V_sc(:,:);


end