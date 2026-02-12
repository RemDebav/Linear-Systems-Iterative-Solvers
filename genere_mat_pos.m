function M=matrice_definie_positive(n)
    N=rand(n);
    M = (N' * N) + eye(n) * 0.1;