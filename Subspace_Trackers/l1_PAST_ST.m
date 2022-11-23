function [W,SEP,eta] = L1_PAST_ST(X_stream,OPTS,U_true)
% Author    : Le Trung Thanh
% School    : University of Orleans, France
% Contact   : thanhle88.tbt@gmail.com


if isfield(OPTS,'lambda'), % Forgetting factor
    lambda = OPTS.lambda;
else lambda = 1;
end
[n, N] = size(X_stream);
r      = size(U_true{1,1},2);


W = eye(n,r);
Z = eye(r);

% mu   = [0  1    0.03    1  0.5    1.0  2.5];
mu = 0.25;

for k = 1 : N
    x     = X_stream(:,k);
    [W,Z] = l1_PAST(x,W,Z,lambda,mu);
    W1    = normalizze(W);
    
    %% Evaluation
    V      = U_true{1,k};
    V      = orth(V);
    SEP(k) = abs(trace(W1'*(eye(n)-V*V')*W1)/trace(W1'*(V*V')*W1));
    eta(k) = sin(subspace(W1,V));
end

end

function [X] = normalizze(X)
%NORMALIZE Normalize the columns (variables) of a data matrix to unit
%Euclidean length.

n = size(X,1);
d = sqrt(sum(X.^2));
d(d == 0) = 1;
X = X./(ones(n,1)*d);

end

function [W P]=l1_PAST(x,W_int,Z_int,B,mu)
P=Z_int;
W=W_int;
y=W'*x;
h=P*y;
g=h./(B+y'*h);

e=x-W*y;
S=(P-g*h');
P=B^-1*Tri(P-g*h');
%%% avant et change P(k-1)
W=W+e*g'+mu*(1-(1/B))*sign(W)*S;
%W=normalizze(W);
%W=orth(W);
W = qrgs(W);
end

function z=Tri(X)
z=triu(X);
z1=z';
z=z+z1-diag(diag(z));
end


function [Q,R] = qrgs(A)
[m,n] = size(A);
% compute QR using Gram-Schmidt
for j = 1:n
    v = A(:,j);
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end
end

