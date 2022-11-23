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
