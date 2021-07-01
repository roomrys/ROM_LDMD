function [rel_energy_Sr, sum_S] = sv_analysis(S)
global showFigs
% plot singular value decay on logarithmic scale
figure('visible', showFigs)
diag_S = diag(S);
semilogy(diag_S ./ S(1, 1), 'LineWidth', 3)
grid on
title('Logarithmic Singular Value Decay Relative to \sigma_1')
ylabel('$log\frac{\sigma_i}{\sigma_1}$', 'interpreter', 'latex')
xlabel('index (i)')
xlim([1 300])
suffix1 = get_suffix1;
saveas(gcf, join([pwd '\Plots\Singular Value Decay\sv_decay' suffix1 '.png'], ''))
% plot relative energy of singular values
sum_S = ones(1, size(S, 1)) * (diag_S);
rel_energy_Sr = cumsum((diag_S)) ./ sum_S;
figure('visible', showFigs)
plot(rel_energy_Sr,'LineWidth',2.0)
title('Relative Energy of Rank Truncated Singular Values')
ylabel({'$\sum_{i = 1:r} \sigma_i / \sum_{j = 1:m} \sigma_j$'},...
    'interpreter', 'latex')
xlabel('Rank Truncation (r)')
xlim([1, 10])
saveas(gcf, join([pwd '\Plots\Singular Value Relative Energy\sv_energy'...
    suffix1 '.png'], ''))
save(join([pwd '\Matrices\Singular Value\S' suffix1 '.mat'], ''), 'diag_S')
