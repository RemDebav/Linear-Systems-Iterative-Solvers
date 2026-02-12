function [X,n]=gauss_implicite(A,b,tol)
    D=diag(diag(A));
    L=tril(A,-1);
    U=triu(A,1);
    N=D+L;
    M=-U;
    B=b;
    X=systeme_triangulaire_inf(N,b);
    n=0;
    while norm(A*X-b)>tol
        X=systeme_triangulaire_inf(N,M*X+B);
        n=n+1;
    end
end