function [e_x, e_y, e_h] = error_analysis(g_3D, Y_pred, AMR, gauges_struct, params)
    global showFigs err rank_trunc
    % find minimum dimensions
    [row_g3D, col_g3D] = size(g_3D);
    [row_Ypred, col_Ypred] = size(Y_pred);
    min_col = min(col_g3D, col_Ypred);
    min_row = min(row_g3D, row_Ypred);
    g_3D = g_3D(1:min_row, 1:min_col);
    Y_pred = Y_pred(1:min_row, 1:min_col);
    
    % total error of each predicted observable
    variables = ['x', 'y', 'h'];
    g_3Dvar = parse_xyh(g_3D);
    clear g_3D
    Y_predvar = parse_xyh(Y_pred);
    clear Y_pred
    if err == 2
        for ii = 1:3
            error.(variables(ii)) = 100*...
                vecnorm(g_3Dvar.(variables(ii)) - Y_predvar.(variables(ii)))...
                ./ vecnorm(g_3Dvar.(variables(ii)));
        end
    else
        for ii = 1:numel(variables)
            g_max.(variables(ii)) = max(g_3Dvar.(variables(ii)), [], 1);
            error.(variables(ii)) = 100*...
                vecnorm(g_3Dvar.(variables(ii)) - Y_predvar.(variables(ii)))...
                ./ vecnorm(g_max.(variables(ii)));
        end
    end
    e_x = error.x; e_y = error.y; e_h = error.h;
    
    % plot error
    g2_fields = fieldnames(AMR.AMR2.gauge_numbers);
    t = gauges_struct.(g2_fields{1})(1:AMR.AMR2.dt_final, 2);
    min_t = min(numel(t), numel(error.x));
    max_y = max([error.x; error.y; error.h], [], 'all');
    min_y = min([error.x; error.y; error.h], [], 'all');
    figure('visible', showFigs)
    line_style = {'b--', 'g-.', 'r:', 'k-.'};
    if abs(max_y - min_y) > 100
        hold on
        for ii = 1:numel(variables)
            semilogy(t(1:min_t), error.(variables(ii))(1:min_t), line_style{ii}, 'LineWidth', 2)
        end
        semilogy([t(params.cutoff_idx), t(params.cutoff_idx)], [0.001, 2*max_y], line_style{ii+1}, 'LineWidth', 2);
        hold off
        grid on
    else
        hold on
        for ii = 1:numel(variables)
            plot(t(1:min_t), error.(variables(ii))(1:min_t), line_style{ii}, 'LineWidth', 2)
        end
        plot([t(params.cutoff_idx), t(params.cutoff_idx)], [0.001, 2*max_y], line_style{ii+1}, 'LineWidth', 2);
        hold off
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
    legend('e_{x}', 'e_{y}', 'e_{h}', 'train data cut-off', 'Location', 'best')
    suffix2 = get_suffix2;
    saveas(gcf, [pwd '\Plots\Error' num2str(err) '\error' suffix2 '.png'])
    save([pwd '\Matrices\Error' num2str(err) '\e_x' suffix2 '.mat'], 'e_x')
    save([pwd '\Matrices\Error' num2str(err) '\e_y' suffix2 '.mat'], 'e_y')
    save([pwd '\Matrices\Error' num2str(err) '\e_eta' suffix2 '.mat'], 'e_h')
    save([pwd '\Matrices\Error' num2str(err) '\t' suffix2 '.mat'], 't')

    % plot largest magnitude observable (if err==3)
    if err ==3
        fig = figure('visible', showFigs);
        for ii = 1:3
            subplot(3, 1, ii)
            plot(t(1:min_t), g_max.(variables(ii))(1:min_t), line_style{ii}, 'LineWidth', 2)
            xlim([0, t(min_t)])
            title(['$\max \vert \xi^n_{' variables(ii) '} \vert$'], 'interpreter', 'latex')
        end
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