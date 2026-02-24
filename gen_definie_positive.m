function A = gen_definie_positive(n, cond_max)
    if nargin < 2, cond_max = 100; end

    R = randn(n);
    A = R' * R;

    lambda_min = min(eig(A));
    if lambda_min < 1e-6
        A = A + (1e-6 - lambda_min + 1) * eye(n);
    end

    if cond(A) > cond_max
        [V, D_eig] = eig(A);
        d_new = linspace(1, cond_max, n)';
        A = V * diag(d_new) * V';
        A = (A + A') / 2;
    end
end
