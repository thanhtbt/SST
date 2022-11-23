function [Ut,rho,eta] = OPAST_Tracking(X,beta,U_tr)
% X:  (n x N) data matrix collecting the N observation (n x 1) vectors 
% beta, 0 < beta<= 1: forgetting factor
% U_tr is the set of true subspaces with time


[n, N] = size(X);
r      = min(size(U_tr{1,1}));

%Initialization
W   = eye(n,r);
Z   = eye(r);
rho = zeros(1,N);
eta = zeros(1,N);

%Processing
for k = 1 : N
    
    y     = W'*X(:,k);
    q     = Z*y/beta;
    GAMA  = 1/(1 + y'*q);
    TAUX  = 1/(norm(q,'fro'))^2*(1/sqrt(1 + (norm(q,'fro'))^2*GAMA^2*(norm(X(:,k),'fro')^2-norm(y,'fro')^2 )) - 1);
    ee    = W*(TAUX*q - GAMA*(1 + TAUX*(norm(q,'fro'))^2)*y) + (1 + TAUX*(norm(q,'fro'))^2)*GAMA*X(:,k);
    Z     = Z/beta - GAMA*(q*q');
    W     = W + ee*q';
    
    % Performance Evaluation
    V       = orth(U_tr{1,k}); % true subspace at time t
    rho(k)  = abs(trace(W'*(eye(n)-V*V')*W)/trace(W'*(V*V')*W));
    eta(k)  = sin(subspace(W,V));
    Ut{1,k} = W;

end