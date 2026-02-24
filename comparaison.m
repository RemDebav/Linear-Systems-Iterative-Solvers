%% Comparaison Jacobi / Gauss-Seidel / SOR
clear; close all; clc;

tailles   = [10, 50, 100, 200, 500, 1000];
tol       = 1e-8;
max_iter  = 5000;
omega_sor = 1.5;

types = {'Triang. sup.', 'Triang. inf.', 'Diag. dominante', ...
         'Tridiagonale', 'Déf. positive', 'Symétrique'};
nb_types = length(types);
nb_n     = length(tailles);

iter_J   = zeros(nb_n, nb_types);
iter_GS  = zeros(nb_n, nb_types);
iter_SOR = zeros(nb_n, nb_types);
time_J   = zeros(nb_n, nb_types);
time_GS  = zeros(nb_n, nb_types);
time_SOR = zeros(nb_n, nb_types);

rng(42);

for t = 1:nb_types
    fprintf('\n=== %s ===\n', types{t});

    for k = 1:nb_n
        n = tailles(k);
        fprintf('  n = %4d ... ', n);

        switch t
            case 1, A = gen_triang_sup(n);
            case 2, A = gen_triang_inf(n);
            case 3, A = gen_diag_dominante(n);
            case 4, A = gen_tridiagonale(n);
            case 5, A = gen_definie_positive(n, 500);
            case 6, A = gen_symetrique(n);
        end

        x_exact = ones(n, 1);
        b = A * x_exact;

        tic; [~, iter_J(k,t)]   = jacobi(A, b, tol, max_iter);       time_J(k,t)   = toc;
        tic; [~, iter_GS(k,t)]  = gauss_seidel(A, b, tol, max_iter); time_GS(k,t)  = toc;
        tic; [~, iter_SOR(k,t)] = sor(A, b, omega_sor, tol, max_iter); time_SOR(k,t) = toc;

        fprintf('J=%4d  GS=%4d  SOR=%4d\n', iter_J(k,t), iter_GS(k,t), iter_SOR(k,t));
    end
end

%% Tracés
couleurs = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]};
noms     = {'Jacobi', 'Gauss-Seidel', sprintf('SOR (\\omega=%.1f)', omega_sor)};

figure('Name', 'Comparaison Jacobi / GS / SOR', 'Position', [50, 100, 1500, 800]);

for t = 1:nb_types
    subplot(2, 3, t);

    yyaxis left;
    semilogy(tailles, iter_J(:,t),   '-o', 'Color', couleurs{1}, 'LineWidth', 1.8, 'MarkerFaceColor', couleurs{1}, 'MarkerSize', 5);
    hold on;
    semilogy(tailles, iter_GS(:,t),  '-s', 'Color', couleurs{2}, 'LineWidth', 1.8, 'MarkerFaceColor', couleurs{2}, 'MarkerSize', 5);
    semilogy(tailles, iter_SOR(:,t), '-^', 'Color', couleurs{3}, 'LineWidth', 1.8, 'MarkerFaceColor', couleurs{3}, 'MarkerSize', 5);
    ylabel('Itérations');
    set(gca, 'YColor', 'k');

    yyaxis right;
    semilogy(tailles, time_J(:,t),   '--o', 'Color', couleurs{1}, 'LineWidth', 1.0, 'MarkerSize', 4);
    semilogy(tailles, time_GS(:,t),  '--s', 'Color', couleurs{2}, 'LineWidth', 1.0, 'MarkerSize', 4);
    semilogy(tailles, time_SOR(:,t), '--^', 'Color', couleurs{3}, 'LineWidth', 1.0, 'MarkerSize', 4);
    ylabel('Temps (s)');
    set(gca, 'YColor', [0.4 0.4 0.4]);
    hold off;

    grid on; xlabel('Taille n'); title(types{t});

    if t == 1
        legend({'J (iter)', 'GS (iter)', 'SOR (iter)', ...
                'J (temps)', 'GS (temps)', 'SOR (temps)'}, ...
               'Location', 'northwest', 'FontSize', 6);
    end
end

sgtitle('Comparaison Jacobi / Gauss-Seidel / SOR', 'FontSize', 14, 'FontWeight', 'bold');
