clc; clear;
%global A a_ij L S sigma tau_ij
load 'armington_data.mat';



options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);

w = fsolve(@fun1,ones(S,1),options);

w = w./w(1);

Y = w.*L;
    
P =  sum(( a_ij.*repmat(A.^(sigma-1),1,S) .* (tau_ij.*repmat(w,1,S)).^(1-sigma) ).^(1./(1-sigma)),1)';
    
X = repmat(A.^(sigma-1),1,S).* repmat(w.^(1-sigma),1,S) .* a_ij .* tau_ij.^(1-sigma) .* repmat( (Y.*P.^(sigma-1))', S,1);
 
