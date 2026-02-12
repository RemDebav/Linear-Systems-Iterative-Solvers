function [X,n]=gauss_explicite(A,b,tol)
    D=diag(diag(A));
    L=tril(A,-1);
    U=triu(A,1);
    invDpL=invtrimatrix(D+L);
    M=-invDpL*U;
    B=invDpL*b;
    X=B;
    n=0;
    while norm(A*X-b)>tol
        X=M*X+B;
        n=n+1;
    end
end