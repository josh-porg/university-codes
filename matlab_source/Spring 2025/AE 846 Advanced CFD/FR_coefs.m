%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Derive the coefficients for a FR-DG scheme
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
k = 3;

format long;
syms x l dl;

% the gauss nodes
if k == 0
    x0(1)=0;
elseif k==1
    x0(1)=-sqrt(1/3);
    x0(2)=-x0(1);
elseif k==2
    x0(1)=-sqrt(3/5);
    x0(2)=0;
    x0(3)=-x0(1);
elseif k == 3
    x0(1)=-sqrt(3/7+2/7*sqrt(6/5));
    x0(2)=-sqrt(3/7-2/7*sqrt(6/5));
    x0(3)=-x0(2);
    x0(4)=-x0(1);
elseif k==4
    x0(1)=-1/3*sqrt(5+2*sqrt(10/7));
    x0(2)=-1/3*sqrt(5-2*sqrt(10/7));
    x0(3)=0;
    x0(4)=-x0(2);
    x0(5)=-x0(1);
else
    disp('not implemeted');
    return;
end

% the Lagrange basis
for i=1:k+1
    l(i)=1;
    for j=1:k+1
        if i ~= j
            l(i) = l(i)*(x-x0(j))/(x0(i)-x0(j));
        end
    end
end

% differential matrix
dl = diff(l);
dmat = zeros(k+1,k+1);
for i=1:k+1
    dmat(i,:)=subs(dl,x,x0(i));
end

% reconstruction coefficients for the values at -1 and 1
recon = zeros(2,k+1);
for i=1:k+1
    recon(1,i)=subs(l(i),x,-1);
    recon(2,i)=subs(l(i),x,1);
end

% Define Legendre polynomial
syms n k

fun = (x^2-1)^n

%P = 1/(2^n*factorial(n))*diff(fun,x,k)

P1 = subs(1/(2^n*factorial(n))*diff(fun, x, 1), n, 1)
P2 = subs(1/(2^n*factorial(n))*diff(fun, x, 2), n, 2)
P3 = subs(1/(2^n*factorial(n))*diff(fun, x, 3), n, 3)
P4 = subs(1/(2^n*factorial(n))*diff(fun, x, 4), n, 4)

rad3 = -1/2*(P3-P2);
drad3 = diff(rad3);
alf_l = double(subs(drad3, x, x0))
alf_r = -flip(alf_l)

