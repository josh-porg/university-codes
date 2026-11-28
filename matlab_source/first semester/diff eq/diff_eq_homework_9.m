syms y(t)
diffEq = diff(diff(y)) + diff(y) -6*y == 12*exp(-2*t);
diffEqHomo = diff(diff(y)) + diff(y) -6*y == 0;
dsolve(diffEq)
dsolve(diffEqHomo)