function [Ut, rho, eta] = FAPI_Tracking(X,beta,U_tr)
% X:  (n x N) data matrix collecting the N observation (n x 1) vectors 
% beta, 0 < beta<= 1: forgetting factor
% U_tr is the set of true subspaces with time


[n, N] = size(X);
r      = size(U_tr{1,1},2);
W      = eye(n,r);
Z      = eye(r);

rho = zeros(1,N); 
eta = zeros(1,N); 

for k = 1:N
    
    y  = W'*X(:,k);
    h  = Z*y;
    g  = h/(beta + y'*h);
    gn = (norm(g,'fro'))^2;

    eps = (norm(X(:,k),'fro'))^2 - (norm(y,'fro'))^2;
    tau = eps/(1 + eps*gn + sqrt(1 + eps*gn) );
    mu  = 1 - tau*gn;
    yp  = mu*y + tau*g;
    hp  = Z'* yp;
    e   = (tau/mu)*(Z*g - (hp'*g)*g);
    Z   = (1/beta)* (Z - g*hp' + e*g');
    ep  = mu* X(:,k) - W*yp;
    W   = W + ep*g';
    
    % Performance Evaluation
    V       = orth(U_tr{1,k}); % true subspace at time t
    rho(k)  = abs(trace(W'*(eye(n)-V*V')*W)/trace(W'*(V*V')*W));
    eta(k)  = sin(subspace(W,V)); 
    Ut{1,k} = W;
end

    
