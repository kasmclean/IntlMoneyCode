clc; clear; % clear data and command window

x0 = [0,0];

tic;
x1 = fsolve(@system1,x0);
t1 = toc;
tic;
x2 = fminunc(@system2,x0); % more brute force than fsolve
t2 = toc;
tic;
xCheck = fminsearch(@system2,x2);
t3 = toc;
tic;
xFrom0 = fminsearch(@system2,x0);
t4 = toc;

function F1 = system1(x)
    F1(1) = exp(-exp(-(x(1)+x(2))))-x(2).*(1+x(1)^2);
    F1(2) = x(1) + x(2) - 1;
end

function F2 = system2(x)
    Y =  exp(-exp(-(x(1)+x(2))))-x(2).*(1+x(1)^2);
    Z = x(1) + x(2) - 1;
    F2 = sum(Y.^2 + Z.^2);
end