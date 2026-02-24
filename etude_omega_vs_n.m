%% Omega optimal en fonction de la taille n
clear; close all; clc;

tailles  = [10, 20, 50, 75, 100, 150, 200, 300, 400, 500,700,1000];
omegas   = 0.05 : 0.025 : 1.975;
tol      = 1e-8;
max_iter = 10000;
cond_max = 500;

nb_n      = length(tailles);
nb_omegas = length(omegas);

omega_opt_exp  = zeros(nb_n, 1);
omega_opt_theo = zeros(nb_n, 1);
iter_at_opt    = zeros(nb_n, 1);
iter_GS        = zeros(nb_n, 1);
iter_Jacobi    = zeros(nb_n, 1);
rho_jacobi     = zeros(nb_n, 1);
iters_map      = NaN(nb_omegas, nb_n);

rng(42);

for k = 1:nb_n
    n = tailles(k);
    fprintf('n = %4d ... ', n);

    A = gen_definie_positive(n, cond_max);
    b = A * ones(n, 1);

    % Rayon spectral de Jacobi + omega théorique (Young)
    D  = diag(diag(A));
    MJ = D \ (D - A);
    rho_jacobi(k)     = max(abs(eig(MJ)));
    omega_opt_theo(k)  = 2 / (1 + sqrt(1 - rho_jacobi(k)^2));

    [~, iter_Jacobi(k)] = jacobi(A, b, tol, max_iter);
    [~, iter_GS(k)]     = gauss_seidel(A, b, tol, max_iter);

    % Balayage de omega
    for w = 1:nb_omegas
        try
            [~, it] = sor(A, b, omegas(w), tol, max_iter);
            if it < max_iter
                iters_map(w, k) = it;
            end
        catch
        end
    end

    [min_it, idx]     = min(iters_map(:, k));
    omega_opt_exp(k)  = omegas(idx);
    iter_at_opt(k)    = min_it;

    fprintf('rho_J=%.4f  w*_theo=%.4f  w*_exp=%.4f  iter*=%d  iterGS=%d  iterJ=%d\n', ...
            rho_jacobi(k), omega_opt_theo(k), omega_opt_exp(k), ...
            iter_at_opt(k), iter_GS(k), iter_Jacobi(k));
end

%% Tracés
figure('Name', 'Omega optimal vs n', 'Position', [100, 150, 1200, 450]);

subplot(1, 2, 1);
plot(tailles, omega_opt_exp,  '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 2, 'MarkerFaceColor', [0.00 0.45 0.74], 'MarkerSize', 7);
hold on;
plot(tailles, omega_opt_theo, '--s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 7);
yline(1, ':', 'Gauss-Seidel', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
hold off;
xlabel('Taille n'); ylabel('\omega^*');
title('\omega^* optimal en fonction de n');
legend({'\omega^* expérimental', '\omega^* théorique (Young)'}, 'Location', 'southeast'); grid on;

subplot(1, 2, 2);
semilogy(tailles, iter_Jacobi, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 2, 'MarkerFaceColor', [0.00 0.45 0.74], 'MarkerSize', 7);
hold on;
semilogy(tailles, iter_GS,     '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 7);
semilogy(tailles, iter_at_opt, '-^', 'Color', [0.47 0.67 0.19], 'LineWidth', 2, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerSize', 7);
hold off;
xlabel('Taille n'); ylabel('Itérations');
title('Nombre d''itérations');
legend({'Jacobi', 'Gauss-Seidel', 'SOR(\omega^*)'}, 'Location', 'northwest'); grid on;

sgtitle(sprintf('Matrice définie positive (cond \\leq %d)', cond_max), 'FontSize', 14, 'FontWeight', 'bold');
