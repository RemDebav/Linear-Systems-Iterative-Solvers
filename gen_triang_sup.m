function U = gen_triang_sup(n, val_min, val_max)
    if nargin < 3, val_max =  10; end
    if nargin < 2, val_min = -10; end

    U = triu(val_min + (val_max - val_min) * rand(n));

    for i = 1:n
        while abs(U(i,i)) < 0.5
            U(i,i) = val_min + (val_max - val_min) * rand();
        end
    end
end
