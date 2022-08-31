function [X_hat,E_hat,L_new,iter] = inexact_alm_modified_rpca_multitwo(D, A, B)

% Nov 2020
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for ajustable Robust PCA.
%
% D - m x n matrix of observations/data (required input)
% InfoMatrix - the embed information matrix.
% lambda - weight on sparse error term in the cost function
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maximum number of iterations
%     - DEFAULT 1000, if omitted or -1.
% Initialize L,E,Y,u, A = InfoMatrix
% while ~converged 
%   minimize (inexactly, update L and E only once)
%     L(L,E,Y,u) = |AL|_* + lambda * |E|_1 + <Y,D-L-E> + mu/2 * |D-L-E|_F^2;
%   Y = Y + \mu * (D - L - E);
%   \mu = \rho * \mu;
% end
%
% Minming Chen, October 2009. Questions? v-minmch@microsoft.com ; 
% Arvind Ganesh (abalasu2@illinois.edu)
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

% addpath PROPACK;

[m,n] = size(D);
lambda = 1/ sqrt(m);
tol = 1e-6;
maxIter = 200;
% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
% [u,s,v]=svd(D,'econ');
% svp = 10;
% L_hat = u(:,1:svp)*s(1:svp,1:svp)*v(:,1:svp)';
X_hat = zeros( m, n);
L_new  = X_hat;
E_hat = zeros( m, n);
mu = 1.25/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.5;          % this one can be tuned,1.15
d_norm = norm(D, 'fro');
iter = 0;
total_svd = 0;
converged = false;
sv = 10;
while ~converged       
    iter = iter + 1;
    % update E_hat
    temp_T = D - L_new + A*X_hat*B -X_hat + (1/mu)*Y;
    E_hat = max(temp_T - lambda/mu, 0);
    E_hat = E_hat+min(temp_T + lambda/mu, 0);
    
    X_hat = D - E_hat + Y/mu;
    if choosvd(n, sv) == 1
        [U,S,V] = lansvd((D - X_hat - E_hat  +Y/mu)*0.5 + A*X_hat*B, sv, 'L');
    else
        [U,S,V] = svd((D - X_hat - E_hat  +Y/mu)*0.5 + A*X_hat*B, 'econ');
    end
%     if choosvd(n, sv) == 1
%         [U,S,V] = lansvd(A*X_hat*B, sv, 'L');
%     else
%         [U,S,V] = svd(A*X_hat*B, 'econ');
%     end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    L_new = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';
    
%     A_hat = slove_generate_nuclear(temp_T,InfoMatrix);
    total_svd = total_svd + 1;
    
    Z = D  - E_hat - X_hat - L_new + A*X_hat*B;
%     Z = D - E_hat -X_hat - L_new - L1 + X_hat*B +A*L1 + (1/mu)*Y;
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    stopCriterion2 = norm(A*X_hat*B-L_new, 'fro') / d_norm;
    if stopCriterion < tol && stopCriterion2 < tol
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
        disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(L_new))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' , stopCriterion1 ' num2str(stopCriterion)...
            ' , stopCriterion2 ' num2str(stopCriterion2)]);
    end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
X_hat = D - E_hat;

