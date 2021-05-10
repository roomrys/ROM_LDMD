function [e_x, e_y, e_eta] = error_analysis(g_3D, Y_pred, AMR2, gauges_struct, params)
    global showFigs err rank_trunc
    % find minimum dimensions
    [row_g3D, col_g3D] = size(g_3D);
    [row_Ypred, col_Ypred] = size(Y_pred);
    min_col = min(col_g3D, col_Ypred);
    min_row = min(row_g3D, row_Ypred);
    % total error of each predicted observable
    if err == 1
        e_x = 100*vecnorm(g_3D(1:3:min_row, 1:min_col) - Y_pred(1:3:min_row, 1:min_col))...
            ./ vecnorm(g_3D(1:min_row, 1:min_col));
        e_y = 100*vecnorm(g_3D(2:3:min_row, 1:min_col) - Y_pred(2:3:min_row, 1:min_col))...
            ./ vecnorm(g_3D(1:min_row, 1:min_col));
        e_eta = 100*vecnorm(g_3D(3:3:min_row, 1:min_col) - Y_pred(3:3:min_row, 1:min_col))...
            ./ vecnorm(g_3D(1:min_row, 1:min_col));
    elseif err == 2
        e_x = 100*vecnorm(g_3D(1:3:min_row, 1:min_col) - Y_pred(1:3:min_row, 1:min_col))...
            ./ vecnorm(g_3D(1:3:min_row, 1:min_col));
        e_y = 100*vecnorm(g_3D(2:3:min_row, 1:min_col) - Y_pred(2:3:min_row, 1:min_col))...
            ./ vecnorm(g_3D(2:3:min_row, 1:min_col));
        e_eta = 100*vecnorm(g_3D(3:3:min_row, 1:min_col) - Y_pred(3:3:min_row, 1:min_col))...
            ./ vecnorm(g_3D(3:3:min_row, 1:min_col));
    else
        g_maxx = max(g_3D(1:3:min_row, 1:min_col), [], 1);
        g_maxy = max(g_3D(2:3:min_row, 1:min_col), [], 1);
        g_maxeta = max(g_3D(3:3:min_row, 1:min_col), [], 1);
        e_x = 100*vecnorm((g_3D(1:3:min_row, 1:min_col) - Y_pred(1:3:min_row, 1:min_col)))...
            ./ g_maxx;
        e_y = 100*vecnorm((g_3D(2:3:min_row, 1:min_col) - Y_pred(2:3:min_row, 1:min_col)))...
            ./ g_maxy;
        e_eta = 100*vecnorm((g_3D(3:3:min_row, 1:min_col) - Y_pred(3:3:min_row, 1:min_col)))...
            ./ g_maxeta;
    end
    g2_fields = fieldnames(AMR2.gauge_numbers);
    t = gauges_struct.(g2_fields{1})(1:AMR2.dt_final, 2);
    min_t = min(numel(t), numel(e_x));
    max_y = max([e_x; e_eta; e_y], [], 'all');
    min_y = min([e_x; e_eta; e_y], [], 'all');
    figure('visible', showFigs)
    if abs(max_y - min_y) > 100
        semilogy(t(1:min_t), e_x(1:min_t), 'b--', ...
        t(1:min_t), e_y(1:min_t), 'g-.', ...
        t(1:min_t), e_eta(1:min_t), 'r:', ...
        [t(params.cutoff_idx), t(params.cutoff_idx)], [0.001, 2*max_y], 'k-.', 'LineWidth', 2);
    grid on
    else
        plot(t(1:min_t), e_x(1:min_t), 'b--', ...
        t(1:min_t), e_y(1:min_t), 'g-.', ...
        t(1:min_t), e_eta(1:min_t), 'r:', ...
        [t(params.cutoff_idx), t(params.cutoff_idx)], [0.001, 2*max_y], 'k-.', 'LineWidth', 2);
    end
    axis([0 t(min_t) min_y max_y])
    ttl = ['\textbf{Relative Error of FOM Wave vs. ROM Wave}'];
    if err == 2
        title({ttl
        ['$e^n = \frac{\Vert \xi^n - \xi_{DMD}^n \Vert_2}{\Vert \xi^n \Vert_2} \qquad r = '...
        num2str(rank_trunc) '$']}, 'interpreter', 'latex')
    else
        title({ttl
        ['$e^n = \frac{\Vert \xi^n - \xi_{DMD}^n \Vert_2}{\max \vert \xi^n \vert} \qquad r = '...
        num2str(rank_trunc) '$']}, 'interpreter', 'latex')
    end
    xlabel('time (s)', 'interpreter', 'latex')
    ylabel('error (\%)', 'interpreter', 'latex')
    legend('e_{x}', 'e_{y}', 'e_{\eta}', 'train data cut-off', 'Location', 'best')
    suffix2 = get_suffix2;
    saveas(gcf, [pwd '\Plots\Error' num2str(err) '\error' suffix2 '.png'])
    save([pwd '\Matrices\Error' num2str(err) '\e_x' suffix2 '.mat'], 'e_x')
    save([pwd '\Matrices\Error' num2str(err) '\e_y' suffix2 '.mat'], 'e_y')
    save([pwd '\Matrices\Error' num2str(err) '\e_eta' suffix2 '.mat'], 'e_eta')
    save([pwd '\Matrices\Error' num2str(err) '\t' suffix2 '.mat'], 't')

    if err ==3
        fig = figure('visible', showFigs);
        subplot(3, 1, 1)
        plot(t(1:min_t), g_maxx(1:min_t), 'b--', 'LineWidth', 2)
        xlim([0, t(min_t)])
        title('$\max \vert \xi^n_{x} \vert$', 'interpreter', 'latex')
        
        subplot(3, 1, 2)
        plot(t(1:min_t), g_maxy(1:min_t), 'g-.', 'LineWidth', 2)
        xlim([0, t(min_t)])
        title('$\max \vert \xi^n_{y} \vert$', 'interpreter', 'latex')
        
        subplot(3, 1, 3)
        plot(t(1:min_t), g_maxeta(1:min_t), 'r:', 'LineWidth', 2)
        xlim([0, t(min_t)])
        title('$\max \vert \xi^n_{\eta} \vert$', 'interpreter', 'latex')
        
        han=axes(fig,'visible','off');
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabel(han,'Largest Magnitude Observable, $\max \vert \xi^n \vert$', 'interpreter', 'latex');
        xlabel(han,'time (s)', 'interpreter', 'latex');
        sgtitle({['\textbf{Largest Magnitude Observable of FOM Wave}']
            ['$\max \vert \xi^n \vert$']}, 'interpreter', 'latex');
        saveas(gcf, [pwd '\Plots\Error' num2str(err) '\max_obsv' suffix2 '.png'])
    end