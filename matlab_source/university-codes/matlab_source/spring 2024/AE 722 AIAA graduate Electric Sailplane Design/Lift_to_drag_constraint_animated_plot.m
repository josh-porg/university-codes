C_D_0 = .0075;
e = .9;
rho = 1.0834; % Density at a+*lt
r = 250;
W_MGTO = 7354.9875; % Maximum Gross Takeoff Weight (N)
b = 20; % wingspan (m)
c = .85; % cord (m)
S = b*c;

A = b^2/S;

L2D = 30
V = convvel(50,"kts","m/s")


W2S = W2S_LiftToDragAtSpeed(L2D, V, C_D_0, A, e, rho)

syms Aspect
fplot(W2S_LiftToDragAtSpeed(L2D, V, C_D_0, Aspect, e, rho), [0,50])
view([90,-90])

L2D = 25
V = convvel(70,"kts","m/s")

W2S = W2S_LiftToDragAtSpeed(L2D, V, C_D_0, A, e, rho)

syms Aspect
fplot(W2S_LiftToDragAtSpeed(L2D, V, C_D_0, Aspect, e, rho), [0,50])
view([90,-90])

L2D = 25
V = convvel(45,"kts","m/s")

W2S = W2S_LiftToDragAtSpeed(L2D, V, C_D_0, A, e, rho)

figure
syms Aspect
fplot(W2S_LiftToDragAtSpeed(L2D, V, C_D_0, Aspect, e, rho), [0,50])
view([90,-90])


figure
drags = linspace(.0025,.0075,100)
for i = 1:length(drags)
    C_D_0 = drags(i)
    syms Aspect
    fplot(W2S_LiftToDragAtSpeed(L2D, V, C_D_0, Aspect, e, rho), [0,50])
    xlim([0,50]); ylim([0,1000]);
    view([90,-90])
    drawnow
    pause(.01)
end
