function [X,n]=jacobi(A,b,tol)
    D=diag(diag(A));
    L=-tril(A,-1);
    U=-triu(A,1);
    invD=diag(1./diag(D));
    M=invD*(L+U);
    B=invD*b;
    X=B;
    n=0;
    while norm(A*X-b)>tol
        X=M*X+B;
        n=n+1;
    end
end