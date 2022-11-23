function [X_stream,U_stream] = online_data_generate(n,T,r,spasity,time_varying,SNR)

U0 = randn(n,r);
I0 = rand(n,r);
I0 = 1 .* (I0 <= 1 - spasity);
U0 = U0 .* I0;
U  = U0;        
noise_level = SNR; %10^(-SNR/20); 
for t = 1 : T
   
    w_t           = randn(r,1);
    x_t           = U * w_t;
    n_t           = randn(n,1);
    X_stream(:,t) = x_t + noise_level * n_t ;
    N             = randn(n,r); 
    epsilon       = time_varying(t);
    if epsilon == 1
         U  = randn(n,r);
         U  = I0 .* U;
    else
         N  = time_varying(t) * N /norm(N,'fro');
         U  = I0 .* (U + N);
    end
   
     U_stream{:,t} = U;
end
end
