function [a b c d] = AtomBalanceLab1(Fuel,ER)
%AE 571 Lab 1

%Atom Balance: Get a,b,c, and d.

if strcmpi(Fuel,'Methane') == 1
    x = 1;
    y = 4;
    astoich = x+y/4;
    a = astoich/ER;
    b = x;
    c = y/2;
    d = a-b-(c/2);
elseif strcmpi(Fuel,'Propane') == 1
    x = 3;
    y = 8;
    astoich = x+y/4;
    a = astoich/ER;
    b = x;
    c = y/2;
    d = a-b-(c/2);
elseif strcmpi (Fuel,'Diesel') == 1
    x = 10.8;
    y = 18.7;
    astoich = x+y/4;
    a = astoich/ER;
    b = x;
    c = y/2;
    d = a-b-(c/2);
elseif strcmpi (Fuel,'Hydrogen') == 1
    x = 0;
    y = 2;
        astoich = x+y/4;
    a = astoich/ER;
    b = x;
    c = y/2;
    d = a-b-(c/2);
end
end