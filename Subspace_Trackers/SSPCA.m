function [U, SEP, eta] = SSPCA(X_stream,OPTS,U_true)
% Code for Streaming Sparse PCA algorithm
% Paper: Yang, Wenzhuo and Huan Xu. "Streaming sparse principal component analysis" 
         ... Proc. ICML, 2015.

% Author    : Le Trung Thanh
% School    : University of Orleans, France
% Contact   : thanhle88.tbt@gmail.com
% Reference : L.T. Thanh et al. "Sparse Subspace Tracking in High Dimensions." 
            ... Proc. IEEE ICASSP, 2022.

[n, N] = size(X_stream);
r      = size(U_true{1,1},2);

if isfield(OPTS,'window'),  % sliding window length
     B = OPTS.window;
else B = round(log(n));
end
if isfield(OPTS,'sparsity'), % sliding window length
    sparsity = OPTS.sparsity;
else sparsity = 0.1;
end

gamma  = round((1-sparsity)*n);


% Evaluation metric
SEP = zeros(1,N);
eta = zeros(1,N);

% Initialization
U   = orth(randn(n,r));

%Processing
k = 1;

while k < N
    if k+B-1 > N
        t_k = k : N;
    else
        t_k = k : k+B-1;
    end
    
    X_k      = X_stream(:,t_k);
    W_k      = U' * X_k;
    S_k_hat  = 1/length(t_k) * X_k * W_k';
    S_k      = row_truncate(S_k_hat,gamma);
    [U,~,~]  = qr(S_k,0);
    
    
    %% Evaluation
    V        =  U_true{1,k};
    V        =  orth(V);
    SEP(t_k) =  abs(trace(U'*(eye(n)-V*V')*U)/trace(U'*(V*V')*U));
    eta(t_k) =  sin(subspace(U,V));
    
    k        =  k + B;
end
end


function S = row_truncate(S_hat,gamma)
    L = size(S_hat,1);
    v = zeros(L,1);
    for ii = 1 : L
        v(ii,1) = norm(S_hat(ii,:));
    end
    [~,idx] = sort(v,'ascend');
    S = S_hat;
    idx_zeros = idx(1:(1-gamma));
    S(idx_zeros,:) = zeros(length(idx_zeros),size(S_hat,2));
end

