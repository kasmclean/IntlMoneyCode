clc; clear; % clear data and command window

load 'armington_data.mat';  % load data for 30 countries

w = ones(S,1); % initial value for w
tolerance = 1; % initial value for tolerance
               % needs to be greater than do loop
update = 0.3; % update is weight of new_w for next iteration's w
              % practical limit of 0.5

% S # of countries
% L_i # labor endowment in country i
% A_i productivity for country i
% tau_ij trade cost between countries i and j
% a_ij exogenous demand shifter
% p_ij price of consuming 1 unit of good from i in country j
% Y_j income of country j
% P_j price index of country j
% X_ij trade between countries i and j
% w_i wages in country i
              
% endogenizing a_ij as f(W)

while tolerance > 0.00001
   
    Y = w.*L; %income of countries in nominal terms
    
    a_ij =ones(S,S)./(abs(repmat(Y,1,S) - repmat(Y',S,1))+ ones(S,S));
    % since a_ij is used to find price level, we should keep calculating 
    % a_ij to nominal values
    
    P =  sum(( a_ij.*repmat(A.^(sigma-1),1,S) .* ...
        (tau_ij.*repmat(w,1,S)).^(1-sigma) ),1)'.^(1./(1-sigma));

    % for technological gap, A_j decreases tau_ij
    
    X = repmat(A.^(sigma-1),1,S).* repmat(w.^(1-sigma),1,S) .* a_ij .* ...
        tau_ij.^(1-sigma) .* repmat( (Y.*P.^(sigma-1))', S,1);
    
    new_w = sum(X,2)./L;
    
    new_w = new_w./new_w(1); % numeraire country is first entry
    
    tolerance = sum(abs(w - new_w)); % could also use sum of squares
    % try to stay close to something nice and boring and well known
    
    disp(['tolerance = ', num2str(tolerance)]); % literally watch tolerance   
    
    w = w.*(1-update) + new_w.*update; % weighted average to update w
    % w = new_w may make system unstable
    
end

