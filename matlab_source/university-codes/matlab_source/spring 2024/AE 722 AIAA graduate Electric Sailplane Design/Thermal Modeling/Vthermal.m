function V_T = Vthermal(r)

V_Tweak = 1.75; % Thermal speed for weak thermal m/s
V_Tstrong = 3.5; % Thermal speed for strong thermal m/s
Grad_nw = 0.025; % Thermal gradient for narrow weak thermals m/s/m (climb speed per radius)
Grad_ns = 0.032; % Thermal gradient for narrow strong thermals m/s/m (climb speed per radius)
Grad_ww = 0.0045; % Thermal gradient for wide weak thermals m/s/m (climb speed per radius)
Grad_ws = 0.006; % Thermal gradient for wide strong thermals m/s/m (climb speed per radius)

V_T = [V_Tweak - Grad_nw * (r - 60); % Thermal type: A1 (narrow, weak) Thermal Speed
       V_Tstrong - Grad_ns * (r - 60); % Thermal type: A2 (narrow, strong) Thermal Speed
       V_Tweak - Grad_ww * (r - 60); % Thermal type: B1 (wide, weak) Thermal Speed
       V_Tstrong - Grad_ws * (r - 60)]; % Thermal type: B2 (wide, strong) Thermal Speed


end