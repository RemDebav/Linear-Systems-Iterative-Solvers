function [X,n]=SOR(A,b,tol,omega)
    D=diag(diag(A));
    L=tril(A,-1);
    U=triu(A,1);
    N=D+omega*L;
    M=-(omega*U+(omega-1)*D);
    B=omega*b;
    X=systeme_triangulaire_inf(N,b);
    n=0;
    while norm(A*X-b)>tol
        X=systeme_triangulaire_inf(N,M*X+B);
        n=n+1;
    end
end