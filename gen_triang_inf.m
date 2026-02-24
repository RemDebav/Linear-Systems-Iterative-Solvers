function L = gen_triang_inf(n, val_min, val_max)
    if nargin < 3, val_max =  10; end
    if nargin < 2, val_min = -10; end

    L = tril(val_min + (val_max - val_min) * rand(n));

    for i = 1:n
        while abs(L(i,i)) < 0.5
            L(i,i) = val_min + (val_max - val_min) * rand();
        end
    end
end
