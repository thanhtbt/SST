function [W,SEP, eta] = SS_FAPI_ST(X_stream,OPTS,U_true)

if isfield(OPTS,'lambda'), % Forgetting factor
     lambda = OPTS.lambda;
else lambda = 1;
end
[n, N] = size(X_stream);
r      = size(U_true{1,1},2);

W = eye(n,r);
Z = eye(r);
Q = eye(r); 
% mu   = [0  1    0.03    1  0.5    1.0  2.5]; 
mu = 0.25;

for k = 1 : N
    x       = X_stream(:,k);
    [W,Z,Q] = SS_FAPI(x, W, Z, Q, lambda, mu,1);
    W1      = normalizze(W*Q);
        
    %% Evaluation
    V      = U_true{1,k};
    V      = orth(V); 
    SEP(k) = abs(trace(W1'*(eye(n)-V*V')*W1)/trace(W1'*(V*V')*W1));
    eta(k) = sin(subspace(W1,V));
end

end

function [W,Z,Q] = SS_FAPI(x, W_init, Z_init,Q_init, Beta,mu,optionmu)
% main algorithm
Z_prev = Z_init;
W_prev = W_init;

y = (W_prev)'*x;
h = Z_prev*y;
g= h./(Beta + y'*h);
e=x-W_prev*y;
tau=(e'*e)/(1+(e'*e)*(g'*g)+sqrt(1+(e'*e)*(g'*g)));
h=Z_prev'*((1-tau*(g'*g))*y+tau*g);
ee=(tau/(1-tau*(g'*g)))*(Z_prev*g-(h'*g)*g);
Z=(1/Beta)*(Z_prev-g*h'+ee*g');
e=(1-tau*(g'*g))*e-tau*W_prev*g;
W = W_prev + e*g';

%Q_init=eye(size(Q_init));

[n, P]=size(W);
M=(W*Q_init);
R=sign(M')*M; %% imposer que Q=Q_init(eye(p)-mu(eps-eps'))
% 
switch optionmu
    case 1
        R=R/sum(sum(R.^2));
        Q=Q_init*(eye(size(W,2))-mu.*R');   
    case 2
        [mu]=0.1/trace(R);
        mu=mu*sum(sum(R.^2));
        R=R/sum(sum(R.^2));
        Q=Q_init*(eye(size(W,2))-mu.*R'); 
    case 3
        R=R/sum(sum(R.^2));
        mu=(sum(abs(M)));%/(trace(R'*R));
         Q=Q_init*(eye(size(W,2))-(ones(P,1)*mu).*R');
    case 4
        R=R/sum(sum(R.^2));
        J= @(x) sum(sum(abs(normalizze(M-x.*M*R'))));
        mu=fminbnd(J,0,100)
        Q=Q_init*(eye(size(W,2))-mu.*R');
    otherwise 
        mu=(sum(abs(M))-1*sqrt(sum(M.^2)))./trace(R'*R);
        Q=Q_init*(eye(size(W,2))-(ones(P,1)*mu).*R');

end
        z=W*Q;
        Q=Q./(ones(P,1)*sqrt(sum(z.^2)));
if any(isnan(Q(:)))||any(isinf(Q(:)))
    error('too big mu')
end

% W=W*Q;
% Q=normalizze(Q);
% Z=(Q'*Z)*inv(Q');
end


function [X] = normalizze(X)
%NORMALIZE Normalize the columns (variables) of a data matrix to unit
%Euclidean length.


n = size(X,1);
d = sqrt(sum(X.^2));
d(d == 0) = 1;
X = X./(ones(n,1)*d);

end


