function A = gen_diag_dominante(n, val_min, val_max)
    if nargin < 3, val_max =  10; end
    if nargin < 2, val_min = -10; end

    A = val_min + (val_max - val_min) * rand(n);

    for i = 1:n
        A(i,i) = 0;
        somme = sum(abs(A(i,:)));
        A(i,i) = somme + 0.5 + (val_max - val_min) * rand();
    end
end
