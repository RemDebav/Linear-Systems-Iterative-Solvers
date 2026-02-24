function [X, nb_iter] = jacobi(A, b, tol, max_iter)
    if nargin < 4, max_iter = 1000; end
    if nargin < 3, tol = 1e-8; end

    D  = diag(diag(A));
    MJ = D \ (D - A);
    cJ = D \ b;
    X  = cJ;

    for k = 1:max_iter
        X_new = MJ * X + cJ;
        if norm(X_new - X) / norm(X_new) < tol
            nb_iter = k;
            return;
        end
        X = X_new;
    end

    nb_iter = max_iter;
    warning('Jacobi : convergence non atteinte en %d itérations.', max_iter);
end
