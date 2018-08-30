clc; clear;

load 'armington_data.mat';

w = ones(S,1); tolerance = 1; update = 0.3;

while tolerance > 0.00001
   
    Y = w.*L;
    
    P =  sum(( a_ij.*repmat(A.^(sigma-1),1,S) .* (tau_ij.*repmat(w,1,S)).^(1-sigma) ),1)'.^(1./(1-sigma));
    
    X = repmat(A.^(sigma-1),1,S).* repmat(w.^(1-sigma),1,S) .* a_ij .* tau_ij.^(1-sigma) .* repmat( (Y.*P.^(sigma-1))', S,1);
    
    new_w = sum(X,2)./L;
    
    new_w = new_w./new_w(1);
    
    tolerance = sum(abs(w - new_w));
    
    disp(['tolerance = ', num2str(tolerance)]);    
    
    w = w.*(1-update) + new_w.*update;
    
end

