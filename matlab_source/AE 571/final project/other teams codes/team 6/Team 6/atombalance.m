function [a, b, c, d, e] = atombalance(x,y,ER)

astoich = x + (y/4);
if ER > 1       %Rich Combustion
elseif ER <= 1   %Lean Combustion
    a = astoich/ER;
    b = x;
    c = (y/2);
    d = a - b - (c/2);
    e = 0;
end
end