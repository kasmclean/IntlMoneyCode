clc; clear; % clear data and command window

load 'armington_data.mat';  % load data for 30 countries

w = ones(S,1); % initial value for w

tolerance = 1; % initial value for tolerance
               % needs to be greater than do loop
update = 0.3; % update is weight of new_w for next iteration's w
              % practical limit of 0.5
tic;
while tolerance > 0.00001
   
    Y = w.*L;
    
    P =  sum(( a_ij.*repmat(A.^(sigma-1),1,S) .* ...
        (tau_ij.*repmat(w,1,S)).^(1-sigma) ),1)'.^(1./(1-sigma));
    
    X = repmat(A.^(sigma-1),1,S).* repmat(w.^(1-sigma),1,S) .* a_ij .* ...
        tau_ij.^(1-sigma) .* repmat( (Y.*P.^(sigma-1))', S,1);
    
    new_w = sum(X,2)./L;
    
    new_w = new_w./new_w(1); % numeraire country is first entry
    
    tolerance = (w-new_w)'*(w-new_w); % could also use sum of squares
    % try to stay close to something nice and boring and well known
    
    disp(['tolerance = ', num2str(tolerance)]); % literally watch tolerance   
    
    w = w.*(1-update) + new_w.*update; % weighted average to update w
    % w = new_w may make system unstable
    
end
t_cmt = toc;

w0 = ones(S,1); % initial value for w
tic;
w1 = fsolve(@armington,w0);
t1 = toc;

tic;
w2 = fminunc(@armington,w0); % more brute force than fsolve
t2 = toc;

tic;
wCheck = fminsearch(@armington,w2);
t3 = toc;

tic;
wFrom0 = fminsearch(@armington,w0);
t4 = toc;

function F1 = armington(w)
    load 'armington_data.mat';
    Y = w.*L;
    load('armington_data.mat', 'sigma');
    P =  sum(( a_ij.*repmat(A.^(sigma-1),1,S) .* ...
        (tau_ij.*repmat(w,1,S)).^(1-sigma) ),1)'.^(1./(1-sigma));
    
    X = repmat(A.^(sigma-1),1,S).* repmat(w.^(1-sigma),1,S) .* a_ij .* ...
        tau_ij.^(1-sigma) .* repmat( (Y.*P.^(sigma-1))', S,1);
    
    new_w = sum(X,2)./L;
    
    new_w = new_w./new_w(1);

    F1 = (w-new_w)'*(w-new_w);  
end


