function F = fun2(w)

global A a_ij L S sigma tau_ij

Y = w.*L;
    
P =  sum(( a_ij.*repmat(A.^(sigma-1),1,S) .* (tau_ij.*repmat(w,1,S)).^(1-sigma) ),1)'.^(1./(1-sigma));
    
X = repmat(A.^(sigma-1),1,S).* repmat(w.^(1-sigma),1,S) .* a_ij .* tau_ij.^(1-sigma) .* repmat( (Y.*P.^(sigma-1))', S,1);
    
wnew = sum(X,2)./L;

F = sum((w - wnew).^2);

end

