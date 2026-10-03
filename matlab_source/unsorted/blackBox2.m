function dX = blackBox2(t, X, Beta)
    x = X(1);
    y = X(2);
    z = X(3);
    
    sigma = Beta(1);
    rho = Beta(2);
    beta = Beta(3);
    
    dX = [
        sigma * (y - x);
        x * (rho - z) - y;
        x * y - beta * z;
    ];
end
