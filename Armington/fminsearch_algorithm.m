clc; clear;

load 'armington_data.mat';

options = optimset('MaxIter',100000,'MaxFunEvals',100000,'Display','iter','PlotFcns', @optimplotfval);

w = fminunc(@fun2,ones(S,1),options);

w = fminsearch(@fun2,w,options);

w = w./w(1);

Y = w.*L;
    
P =  sum(( a_ij.*repmat(A.^(sigma-1),1,S) .* (tau_ij.*repmat(w,1,S)).^(1-sigma) ).^(1./(1-sigma)),1)';
    
X = repmat(A.^(sigma-1),1,S).* repmat(w.^(1-sigma),1,S) .* a_ij .* tau_ij.^(1-sigma) .* repmat( (Y.*P.^(1-sigma))', S,1);

