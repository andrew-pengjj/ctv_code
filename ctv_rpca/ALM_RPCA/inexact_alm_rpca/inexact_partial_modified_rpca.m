function [R_hat,E_hat,iter] = inexact_partial_modified_rpca(D, G)

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
%     L(L,E,Y,u) = |(I-GG^*)L|_* + lambda * |E|_1 + <Y,D-L-E> + mu/2 * |D-L-E|_F^2;
%           Y = Y + \mu * (D - L - E);
%           \mu = \rho * \mu;
% end
%
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

% addpath PROPACK;

[m,n] = size(D);
lambda = 1 / sqrt(max(m,n));
tol = 1e-6;
maxIter = 100;
% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

E_hat = zeros( m, n);
% [u,s,v] = svd((eye(m)-G*G')*D, 'econ');
% L_hat   = u(:,1:10)*s(1:10,1:10)*v(:,1:10)';
L_hat = zeros( m, n);
X_hat = L_hat'*G; 

mu = 1.25/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.25;          % this one can be tuned,1.15
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
sv = 10;
while ~converged       
    iter = iter + 1;
    % update E_hat
    temp_T = D - G*X_hat' - L_hat + (1/mu)*Y;
    E_hat  = max(temp_T - lambda/mu, 0);
    E_hat  = E_hat+min(temp_T + lambda/mu, 0);
    
    % update L_hat
    if choosvd(n, sv) == 1
        [U,S,V] = lansvd((eye(m)-G*G')*(D - E_hat + (1/mu)*Y), sv, 'L');
    else
        [U,S,V] = svd((eye(m)-G*G')*(D - E_hat + (1/mu)*Y), 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    L_hat = (U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)');
    
    % update X_hat
    X_hat = (D - E_hat+(1/mu)*Y)'*G;
    
    total_svd = total_svd + 1;
    Z = D - L_hat - E_hat - G*X_hat';
   
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol && iter >20
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
        disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(L_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
    end    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
    
end
R_hat = D-E_hat;
