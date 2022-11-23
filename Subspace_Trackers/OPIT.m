function [U, SEP, rho] = OPIT(X_stream,OPTS,U_true)
% Author    : Le Trung Thanh
% School    : University of Orleans, France
% Contact   : thanhle88.tbt@gmail.com
% Reference : L.T. Thanh et al. "Sparse Subspace Tracking in High Dimensions." 
            ... Proc. IEEE ICASSP, 2022.

[n, N] = size(X_stream);
r      = size(U_true{1,1},2);

if isfield(OPTS,'lambda'), % forgetting factor
     lambda = OPTS.lambda;
else lambda = 0.99; 
end
if isfield(OPTS,'window'), % sliding window length
     B = OPTS.window;
else
    if n < 1e3 B = 1; 
    else  B = round(2*log(n));
    end
end
if isfield(OPTS,'method'), % normalization or orthonormalization?
     method  = OPTS.method;
else
    if n < 1e3 method = 'normalization';
    else  method = 'orthonormalization';
    end
end

if isfield(OPTS,'omega'), % sparsity
     omega = OPTS.omega;
     m     = round((1-omega)*n);
else m     = round(10*r*log(n));
end

% Evaluation metric
SEP = zeros(1,N);
rho = zeros(1,N);

% Initialization
U0  = orth(randn(n,r));
U   = U0;
S_k = zeros(n,r);
E_k = eye(r);


%% Processing
k = 1;
while k < N
    if k+B-1 > N
        t_k = k:N;
    else
        t_k = k:k+B-1;
    end
    X_k     = X_stream(:,t_k);
    U_old   = U;
    Z_k     = U_old' * X_k;
    S_old   = S_k;
    S_k     = lambda * (B*(k-1))/(B*k) * S_k * E_k + 1/(B*k) * X_k * Z_k';
    Idx_max = maxk(S_k,m,'ComparisonMethod','abs');
    for r_i = 1 : r
        Skx           = S_k(:,r_i);
        idx_max       = Idx_max(:,r_i);
        Rkx(~idx_max) = 0;
        S_k(:,r_i)    = Skx;
    end
    if strcmp(method,'normalization')
        U = S_k / norm(S_k);  % normalization
    else
        [U,~,~] = qr(S_k,0);  % QR-orthonormalization
    end
    E_k = U_old' * U;
    
    %% Evaluation
    V        =  U_true{1,k};
    V        =  orth(V);
    SEP(t_k) =  abs(trace(U'*(eye(n)-V*V')*U)/trace(U'*(V*V')*U));
    rho(t_k) =  sin(subspace(U,V));
    k        =  k + B;
    
end
SEP(1) =  abs(trace(U0'*(eye(n)-V*V')*U0)/trace(U0'*(V*V')*U0));
rho(1) =  sin(subspace(U0,V));
end
