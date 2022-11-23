function [U, SEP, eta] = Oja(X_stream,U_true)
% Author    : Le Trung Thanh
% School    : University of Orleans, France
% Contact   : thanhle88.tbt@gmail.com
% Reference : L.T. Thanh et al. "Sparse Subspace Tracking in High Dimensions." 
            ... Proc. IEEE ICASSP, 2022.

[n, N] = size(X_stream);
r      = size(U_true{1,1},2);

SEP = zeros(1,N);
eta = zeros(1,N);

%Initialization

U = orth(randn(n,r));

mu = 1;

for t = 1 : N
    beta  = 1/t; 
    x_t   = X_stream(:,t);
    z_t   = U'*x_t;
    S_t   = U + beta * mu * x_t * z_t';
    [U,~,~] = qr(S_t,0);
    
    
    %% Evaluation
    V      =  U_true{1,t};
    V      =  orth(V);
    SEP(t) =  abs(trace(U'*(eye(n)-V*V')*U)/trace(U'*(V*V')*U));
    eta(t) =  sin(subspace(U,V));
   
end


