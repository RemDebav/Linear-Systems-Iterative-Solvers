%% Influence de omega sur SOR
clear; close all; clc;

n         = 100;
omegas    = 0.05 : 0.05 : 1.95;
tol       = 1e-8;
max_iter  = 10000;

types = {'Triang. sup.', 'Triang. inf.', 'Diag. dominante', ...
         'Tridiagonale', 'Déf. positive', 'Symétrique'};
nb_types  = length(types);
nb_omegas = length(omegas);

iters  = NaN(nb_omegas, nb_types);
temps  = NaN(nb_omegas, nb_types);
rayons = NaN(nb_omegas, nb_types);

%% Génération des matrices
rng(42);
matrices = cell(1, nb_types);
seconds  = cell(1, nb_types);

for t = 1:nb_types
    switch t
        case 1, A = gen_triang_sup(n);
        case 2, A = gen_triang_inf(n);
        case 3, A = gen_diag_dominante(n);
        case 4, A = gen_tridiagonale(n);
        case 5, A = gen_definie_positive(n, 500);
        case 6, A = gen_symetrique(n);
    end
    matrices{t} = A;
    seconds{t}  = A * ones(n, 1);
end

%% Balayage de omega
for t = 1:nb_types
    A = matrices{t};
    b = seconds{t};

    D = diag(diag(A));
    L = -(tril(A, -1));
    U = -(triu(A, 1));

    fprintf('=== %s ===\n', types{t});

    for w = 1:nb_omegas
        omega = omegas(w);

        DwL   = D - omega * L;
        M_sor = DwL \ ((1 - omega) * D + omega * U);
        rayons(w, t) = max(abs(eig(M_sor)));

        try
            tic;
            [~, it] = sor(A, b, omega, tol, max_iter);
            temps(w, t) = toc;

            if it < max_iter
                iters(w, t) = it;
            else
                iters(w, t) = NaN;
                temps(w, t) = NaN;
            end
        catch
            iters(w, t) = NaN;
            temps(w, t) = NaN;
        end
    end
    fprintf('  terminé.\n');
end

%% Figure 1 : par type de matrice
couleurs_type = lines(nb_types);

figure('Name', 'Influence de omega par type', 'Position', [30, 80, 1600, 750]);

for t = 1:nb_types
    subplot(2, 3, t);

    yyaxis left;
    plot(omegas, iters(:, t), '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', [0.00 0.45 0.74]);
    ylabel('Itérations');
    set(gca, 'YColor', [0.00 0.45 0.74]);

    yyaxis right;
    plot(omegas, rayons(:, t), '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', [0.85 0.33 0.10]);
    hold on;
    yline(1, '--', '\rho = 1', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    hold off;
    ylabel('\rho(M_\omega)');
    set(gca, 'YColor', [0.85 0.33 0.10]);

    xlabel('\omega'); title(types{t}); grid on;

    if t == 1
        legend({'Itérations', '\rho(M_\omega)', 'Seuil \rho=1'}, 'Location', 'north', 'FontSize', 7);
    end

    % Omega optimal
    [min_it, idx] = min(iters(:, t));
    if ~isnan(min_it)
        yyaxis left; hold on;
        plot(omegas(idx), min_it, 'p', 'Color', [0.47 0.67 0.19], 'MarkerSize', 14, 'MarkerFaceColor', [0.47 0.67 0.19]);
        text(omegas(idx), min_it * 1.3, sprintf('\\omega^*=%.2f', omegas(idx)), 'FontSize', 7, 'Color', [0.47 0.67 0.19], 'FontWeight', 'bold');
        hold off;
    end
end

sgtitle(sprintf('Influence de \\omega sur SOR  (n = %d)', n), 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 2 : superposition tous types
figure('Name', 'Comparaison tous types', 'Position', [30, 80, 1400, 500]);

subplot(1, 3, 1);
for t = 1:nb_types
    semilogy(omegas, iters(:, t), '-', 'Color', couleurs_type(t,:), 'LineWidth', 1.5); hold on;
end
hold off;
xlabel('\omega'); ylabel('Itérations'); title('Itérations');
legend(types, 'Location', 'northwest', 'FontSize', 7); grid on;

subplot(1, 3, 2);
for t = 1:nb_types
    semilogy(omegas, temps(:, t), '-', 'Color', couleurs_type(t,:), 'LineWidth', 1.5); hold on;
end
hold off;
xlabel('\omega'); ylabel('Temps (s)'); title('Temps');
legend(types, 'Location', 'northwest', 'FontSize', 7); grid on;

subplot(1, 3, 3);
for t = 1:nb_types
    plot(omegas, rayons(:, t), '-', 'Color', couleurs_type(t,:), 'LineWidth', 1.5); hold on;
end
yline(1, '--k', '\rho = 1', 'LineWidth', 1.2);
hold off;
xlabel('\omega'); ylabel('\rho(M_\omega)'); title('Rayon spectral');
legend([types, {'Seuil \rho=1'}], 'Location', 'northwest', 'FontSize', 7); grid on;

sgtitle(sprintf('Influence de \\omega — tous types  (n = %d)', n), 'FontSize', 14, 'FontWeight', 'bold');
