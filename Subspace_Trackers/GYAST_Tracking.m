function [Ut,rho,eta] = GYAST_Tracking(X,Beta,U_tr)
% Author    : Le Trung Thanh
% School    : University of Orleans, France
% Contact   : thanhle88.tbt@gmail.com

% X:  (n x N) data matrix collecting the N observation (n x 1) vectors
% beta, 0 < beta<= 1: forgetting factor
% U_tr is the set of true subspaces with time

[n, N] = size(X);
p      = min(size(U_tr{1,1}));

rho = zeros(1,N);
eta = zeros(1,N);

sigma = 0.0001;
opts.LT=true;
opts1.UT=true;

Cy2 = sigma * eye(p,'double');
Um2 = eye(n,p,'double');
Cyprim2 = zeros(p+1,p+1,'double');
qn2 = eye(p+1,1);

    for k=1:N
        % Input vector
        x=X(:,k);
        
        if k==1
            varN2=(norm(x))^2/n;
        end
        y=Um2'*x;
        Sigmma = sqrt(x'*x-y'*y);
        Cprimyy = Beta*Cy2+y*y';
        Cyprim2(1:p,1:p) = Cprimyy;
        z = Sigmma*y;
        Cyprim2(1:p,p+1) = z;
        Cyprim2(p+1,1:p) = z';
        GAMMA = Beta*varN2+Sigmma^2;
        Cyprim2(p+1,p+1) = GAMMA;
        Cyprim2=(Cyprim2+Cyprim2')/2;
        UT=chol(Cyprim2);
        for gg=1:6
  % %            UT'*bbb=qn;
                bbb=linsolve(UT',qn2,opts);
  % %            UT*aaa=bbb;
                aaa=linsolve(UT,bbb,opts1);
                qn2=aaa/norm(aaa);
        end
        %%%%%
%         [qn2,~]=eigs(Cyprim2,1,'sm');
        q = qn2(1:p);
        r = qn2(p+1);
        l = 1 / (1+abs(r));
        K = r / (Sigmma*abs(r));
        eprime = Um2*(l*q-K*y)+K*x;
        Um2 = Um2-eprime*q';
        qprim = Cprimyy * q;
        C = (r*z/abs(r)+l*qprim)*q';
        c = l^2*q'*qprim + GAMMA + 2*l*real(r*q'*z)/abs(r);
        Cy2 = Cprimyy + c*q*q' - C - C'; 
        Cy2 = (Cy2+Cy2')/2;
%         Cy2 = Beta*Cy2 + (y*y');
        LAMBDAm = trace(Cyprim2)-trace(Cy2);
        varN2 = min([LAMBDAm varN2]);
        
        
        
        V      = orth(U_tr{1,k}); % true subspace at time t
        rho(k)  = abs(trace(Um2'*(eye(n)-V*V')*Um2)/trace(Um2'*(V*V')*Um2));
        eta(k)  = sin(subspace(Um2,V)); 
        Ut{1,k} = Um2;
    
    
        
    end
   
end