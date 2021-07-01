function err_struct = error_analysis(G, Y_pred, chosen, AMR, gauges_struct, params)
    global showFigs err rank_trunc
    % find minimum dimensions
    [row_g3D, col_g3D] = size(G);
    [row_Ypred, col_Ypred] = size(Y_pred);
    min_col = min(col_g3D, col_Ypred);
    min_row = min(row_g3D, row_Ypred);
    G = G(1:min_row, 1:min_col);
    Y_pred = Y_pred(1:min_row, 1:min_col);
    
    % total error of each predicted observable
    % TO2D: user selected variables only here!
    Gvar = parse_xyh(G, chosen);
    clear G
    Y_predvar = parse_xyh(Y_pred, chosen);
    clear Y_pred
    num_chosen = numel(chosen);
    if err == 2
        for ii = 1:num_chosen
            err_struct.(chosen(ii)) = 100*...
                vecnorm(Gvar.(chosen(ii)) - Y_predvar.(chosen(ii)))...
                ./ vecnorm(Gvar.(chosen(ii)));
        end
    else
        for ii = 1:num_chosen
            g_max.(chosen(ii)) = max(Gvar.(chosen(ii)), [], 1);
            err_struct.(chosen(ii)) = 100*...
                vecnorm(Gvar.(chosen(ii)) - Y_predvar.(chosen(ii)))...
                ./ vecnorm(g_max.(chosen(ii)));
        end
    end
    suffix2 = get_suffix2;
    save(join([pwd '\Matrices\Error' num2str(err) '\err_struct' suffix2 '.mat'], ''), 'err_struct');
    
    % plot error
    g2_fields = fieldnames(AMR.AMR2.gauge_numbers);
    t = gauges_struct.(g2_fields{1})(1:AMR.AMR2.dt_final, 2);
    min_t = min(numel(t), numel(err_struct.(chosen(1))));
    error_names = strcat(repmat("e_{", 1, num_chosen), chosen, "}");
    max_y = -Inf; min_y = Inf;
    for ii = 1:num_chosen
        max_y = max([max_y; max(err_struct.(chosen(ii)))], [], 'all');
        min_y = min([min_y; min(err_struct.(chosen(ii)))], [], 'all');
    end
        
%     max_y = max([err_struct.x; err_struct.y; err_struct.h], [], 'all');
%     min_y = min([err_struct.x; err_struct.y; err_struct.h], [], 'all');
    figure('visible', showFigs)
    line_style = {'b--', 'g-.', 'r:', 'k-.'};
    hold on
    if abs(max_y - min_y) > 100
        for ii = 1:num_chosen
            semilogy(t(1:min_t), err_struct.(chosen(ii))(1:min_t), line_style{ii}, 'LineWidth', 2)
        end
        semilogy([t(params.cutoff_idx), t(params.cutoff_idx)], [0.001, 2*max_y], line_style{ii+1}, 'LineWidth', 2);
        grid on
    else
        for ii = 1:num_chosen
            plot(t(1:min_t), err_struct.(chosen(ii))(1:min_t), line_style{ii}, 'LineWidth', 2)
        end
        plot([t(params.cutoff_idx), t(params.cutoff_idx)], [0.001, 2*max_y], line_style{end}, 'LineWidth', 2);
    end
    legend([error_names, "train data cut-off"], 'Location', 'best')
    hold off
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
    saveas(gcf, join([pwd '\Plots\Error' num2str(err) '\error' suffix2 '.png'], ''))
    save(join([pwd '\Matrices\Error' num2str(err) '\t' suffix2 '.mat'], ''), 't')

    % plot largest magnitude observable (if err==3)
    if err ==3
        fig = figure('visible', showFigs);
        for ii = 1:3
            subplot(3, 1, ii)
            plot(t(1:min_t), g_max.(chosen(ii))(1:min_t), line_style{ii}, 'LineWidth', 2)
            xlim([0, t(min_t)])
            title(['$\max \vert \xi^n_{' chosen(ii) '} \vert$'], 'interpreter', 'latex')
        end
        han=axes(fig,'visible','off');
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabel(han,'Largest Magnitude Observable, $\max \vert \xi^n \vert$', 'interpreter', 'latex');
        xlabel(han,'time (s)', 'interpreter', 'latex');
        sgtitle({['\textbf{Largest Magnitude Observable of FOM Wave}']
            ['$\max \vert \xi^n \vert$']}, 'interpreter', 'latex');
        saveas(gcf, join([pwd '\Plots\Error' num2str(err) '\max_obsv' suffix2 '.png'], ''))
    end