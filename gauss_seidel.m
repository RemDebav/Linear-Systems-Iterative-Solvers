function [X, nb_iter] = gauss_seidel(A, b, tol, max_iter)
    if nargin < 4, max_iter = 1000; end
    if nargin < 3, tol = 1e-8; end

    DL  = tril(A);
    U   = -(triu(A, 1));
    MGS = DL \ U;
    cGS = DL \ b;
    X   = cGS;

    for k = 1:max_iter
        X_new = MGS * X + cGS;
        if norm(X_new - X) / norm(X_new) < tol
            nb_iter = k;
            return;
        end
        X = X_new;
    end

    nb_iter = max_iter;
    warning('Gauss-Seidel : convergence non atteinte en %d itérations.', max_iter);
end
