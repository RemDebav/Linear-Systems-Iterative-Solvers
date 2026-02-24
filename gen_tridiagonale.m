function A = gen_tridiagonale(n, val_min, val_max, dominante)
    if nargin < 4, dominante = true; end
    if nargin < 3, val_max   =  10;  end
    if nargin < 2, val_min   = -10;  end

    sub  = val_min + (val_max - val_min) * rand(n-1, 1);
    sup  = val_min + (val_max - val_min) * rand(n-1, 1);
    diag_vals = val_min + (val_max - val_min) * rand(n, 1);

    A = diag(diag_vals) + diag(sub, -1) + diag(sup, 1);

    if dominante
        for i = 1:n
            somme = sum(abs(A(i,:))) - abs(A(i,i));
            A(i,i) = sign(A(i,i)) * (somme + 0.5 + abs(val_max) * rand());
            if abs(A(i,i)) < 0.5
                A(i,i) = somme + 1;
            end
        end
    end
end
