function [u, p, E, rho, c] = getFlowProperties(Q, gamma) % checked - not confident (ok?)
    % get thermodynamic properties of flow from state vector Q
    rho = Q(1,:);
    u = Q(2,:) ./ Q(1,:);
    p = (Q(3,:) - 1/2 .* Q(2,:).^2./Q(1,:)) * (gamma - 1); % currently p is going negative sometimes
    E = Q(3,:);
    c = sqrt(gamma*p./rho); % given by josh sutton and jackson (this is correct)
    
    % FOR DEBUG
    if sum(rho < 0) > 0
        disp("negative rho There is a bug")
    end
    if sum(p < 0) > 0
        disp("negative p There is a bug")
    end
    if sum((gamma*p./rho) < 0) > 0
        disp("complex spd sound There is a bug")
    end
    %c = 1; % FOR DEBUG
end