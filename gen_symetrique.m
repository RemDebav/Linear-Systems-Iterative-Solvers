function A = gen_symetrique(n, val_min, val_max)
    if nargin < 3, val_max =  10; end
    if nargin < 2, val_min = -10; end

    R = val_min + (val_max - val_min) * rand(n);
    R = R - diag(diag(R));
    A = (R + R') / 2;

    for i = 1:n
        A(i,i) = sum(abs(A(i,:))) + 0.5 + abs(val_max) * rand();
    end
end
