function plot_EV(D)
global showFigs
figure('visible', showFigs)
plot(real(D(~(abs(D) > 1))), imag(D(~(abs(D) > 1))), 'bo', 'LineWidth', 1)
hold on
plot(real(D(abs(D) > 1)), imag(D(abs(D) > 1)), 'ro', 'LineWidth', 1)
title({['Eigenvalues of $K$ where $\xi_{n+1} = K \xi_n$']
    ['Number $|D| \leq 1$: ' num2str(sum(abs(diag(D)) <= 1)) '; Number $|D| > 1$: ' num2str(numel(D(abs(D) > 1)))]}, 'interpreter', 'latex')
xlabel('Re\{D\}', 'interpreter', 'latex')
ylabel('Im\{D\}', 'interpreter', 'latex')
theta = 0:0.01:(2*pi);
plot(cos(theta), sin(theta), 'k--')
hold off
suffix2 = get_suffix2;
legend('$|D| \leq 1$', '$|D| > 1$', 'interpreter', 'latex', 'Location', 'best')
saveas(gcf, [pwd '\Plots\LDMD Spectral Analysis\ev_uc_K' suffix2 '.png'])
xlim([min(min(real(D)))*1.10, max(max(real(D)))*1.10])
ylim([min(min(imag(D)))*1.10, max(max(imag(D)))*1.10])
legend('$|D| \leq 1$', '$|D| > 1$', 'interpreter', 'latex', 'Location', 'best')
saveas(gcf, [pwd '\Plots\LDMD Spectral Analysis\ev_K' suffix2 '.png'])
