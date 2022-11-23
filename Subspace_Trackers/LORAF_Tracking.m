function [Ut,rho,eta] = LORAF_Tracking(X,beta,U_tr)
% X:  (n x N) data matrix collecting the N observation (n x 1) vectors
% beta, 0 < beta<= 1: forgetting factor
% U_tr is the set of true subspaces with time

[n, N] = size(X);
r      = min(size(U_tr{1,1}));

%Initialization
D   = eye(n,r);
%D   = orth(randn(n,r));
%Z   = eye(n);
%P   = eye(r);

S = 0.00001*eye(r);
v = zeros(r,1);

% performance assessment
rho = zeros(1,N);
eta = zeros(1,N);
%Processing
% for k = 1:N
%     x = X(:,k);
%     y = D'*x;
%     % caluclate weight
%     x2 = x'*x;
%     y2 = y'*y;
% %     trCx = beta*trCx + x2;
% %     trCy = beta*trCy + y2;
% %     sig2n = (trCx-trCy)/(n-r);
% %     ytild = Siv*y;
%     ZZ = x2-y2;
% %     w = 1/(ytild'*y+ ZZ/sig2n);
% %    w = l*(1/(ytild'*y+ ZZ/sig2n)) + (1-l)*(1/x2);
% 
% %  % caluclate weight
% %     uu = Z*x;
% %     delt = x'*uu;
% %     vv = uu/(beta + delt);
% %     Z = (Z-uu*vv')/beta;
% %     w = 1/(delt); wstore(k) = w;
%    
%    % w = 1;
%     
%     % calculate subspace 
%     etav = S*v;
%     Xv = beta*S-2*beta*etav*v' + (y*y');
%     b = pinv(Xv')*sqrt(ZZ)*y;
%     vt = 4*(b'*b +1);
%     ph2 = 1/2 + 1/sqrt(vt);
%     gam = (1-2*ph2)/(2*sqrt(ph2));
%     xi = sqrt(ph2)/sqrt(ZZ);
%     v = gam*b;
%     S = Xv - (1/xi)*(v*y');
%     
%     %
%     uv = xi*y - v;
%     e = xi*x - D*uv;
%     D = D-2*e*v';
%     
%     rho(k)=  abs(trace(D'*(eye(n)-V*V')*D)/trace(D'*(V*V')*D));
%     eta(k)=  norm(D' * D - eye(r))^2;
% end

for k = 1:N
    z = X(:,k);
    h = D'*z;
    Z = z'*z - h'*h;
    if Z == 0
        print('Oh la la');
    end
    u = S*v;
%    XX = beta*S + 2*beta*u*v' + h*h';
    XX = beta*S - 2*beta*u*v' + h*h';
    b = sqrt(Z)*pinv(XX')*h;
    bet = 4*(b'*b + 1);
    p2 = 1/2 + 1/sqrt(bet);
    gam = (1-2*p2)/(2*sqrt(p2));
    del = sqrt(p2)/sqrt(Z);
    v = gam*b;
    S = XX-(1/del)*v*h';
    w = del*h-v;
    e = del*z - D*w;
    D = D - 2*e*v';
    
    % performance 
    V       = orth(U_tr{1,k}); % true subspace at time t
    rho(k)  = abs(trace(D'*(eye(n)-V*V')*D)/trace(D'*(V*V')*D));
    eta(k)  = sin(subspace(D,V)); 
    Ut{1,k} = D;
end

