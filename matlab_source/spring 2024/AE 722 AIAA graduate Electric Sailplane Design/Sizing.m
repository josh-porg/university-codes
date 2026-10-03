% Sizing to TO
clc;
close all
W_to_max=2180;        %Max takeoff weight in lbs from intial sizing
We=1020;                %empty weight in lbs
S=174;                %Average wing area from STAMPED in ft^2
S_to_max=1640;       %Max takeoff distance in meters from CS 22 (= 500 m)
S_tog_max = S_to_max/1.6;    %From Fig 3.4 pg 96, S_tog from S_to
TOP23=86;                 %Eq 3.6, in lbs^2/(ft^2hp)
CL_max_to= linspace(1.4,2,3);
sigma=1;                     %atmospheric denisty ratio (1 at sea level)
WS_to=linspace(We/S,W_to_max*1.5/S,500);
WP_to1=TOP23*sigma*CL_max_to(:,1)./(WS_to);
WP_to2=TOP23*sigma*CL_max_to(:,2)./(WS_to);
WP_to3=TOP23*sigma*CL_max_to(:,3)./(WS_to);
plot(WS_to,WP_to1,WS_to,WP_to2,WS_to,WP_to3)
hold on
xline(W_to_max/S)
grid on
ylabel('Power Loading (lbs/hp)'), xlabel('Wing Loading (lbs/ft^2)')
legend('C_L_T_O=1.4','C_L_T_O=1.7','C_L_T_O=2.0','Wing Loading at MTOW')
%% Climb sizing
% Drag Polars
clc
Swet=480;   % Estimate from fig 3.22
cf=.005;
AR=22;       % DG1000
e=.9;
CL_max=linspace(1.2,1.8,3);
CL=linspace(0,3,50);
f_graph=Swet*(CL/46.5-(CL/(pi*AR*e)));
f=Swet*(CL_max/46.5-(CL_max/(pi*AR*e)));
Cd_clean=f_graph/Swet+CL.^2/(pi*AR*.9);    %clean
Cd_fl=.01+f_graph/Swet+CL.^2/(pi*AR*.85);       %Takeoff flaps
Cd_flla=.065+f_graph/Swet+CL.^2/(pi*AR*.8);    %Takeoff flaps gear down
plot(CL,Cd_clean,CL,Cd_fl,CL,Cd_flla)
xlabel('Coefficient of Lift, CL (~)'), ylabel('Coefficient of Drag, CD(~)')
xline(CL_max(:,1),'b--'), xline(CL_max(:,2),'k--'), xline(CL_max(:,3),'r--')
grid on
legend('Clean Config.','TO Flaps','TO Flaps & Landing Gear','C_L_m_a_x=1.2','C_L_m_a_x=1.5','C_L_m_a_x=1.8')

%% Climb  (dh/dt)
clc;
rho=.0023; 
sigma=.9761;
WS_climb=linspace(We/S,W_to_max/S,10);
eta_p=.8;
RCP=[295 600 1200]/33000;
Cdo=f(:,2)/Swet+.01;    %clean plus flaps
Vclimb=sqrt(2*WS_climb/(rho*CL_max_to(:,3)));      %CLmax with flaps in best case 
CL_CD_32max=1.345*(AR*e)^.75/(Cdo^.25);
RCP_max=eta_p/WP_to3(:,250)-sqrt(W_to_max/S)/(19*sqrt(sigma)*CL_CD_32max);
W_P_parameters=eta_p./(sqrt(W_to_max/S)./(19*sqrt(sigma)*CL_CD_32max)+RCP);

