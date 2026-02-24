function [X, nb_iter] = sor(A, b, omega, tol, max_iter)
    if nargin < 5, max_iter = 1000; end
    if nargin < 4, tol = 1e-8; end
    if nargin < 3, omega = 1.5; end

    D = diag(diag(A));
    L = -(tril(A, -1));
    U = -(triu(A, 1));

    DwL  = D - omega * L;
    Msor = DwL \ ((1 - omega) * D + omega * U);
    cSOR = omega * (DwL \ b);
    X    = cSOR;

    for k = 1:max_iter
        X_new = Msor * X + cSOR;
        if norm(X_new - X) / norm(X_new) < tol
            nb_iter = k;
            return;
        end
        X = X_new;
    end

    nb_iter = max_iter;
    warning('SOR : convergence non atteinte en %d itérations.', max_iter);
end
