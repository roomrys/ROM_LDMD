% Physics Aware DMD
% note: each field in g_struct has following format--
% # level, time, q[  1  2  3], eta, aux[]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%TO2D: add variable for which observables to use (all 3 or only 2 -- which 2?)
global sim_number showFigs err isScaled all_amr art_stab chosen_v
chosen_v = ["x", "y", "h"];
chosen_v = chosenObservables(chosen_v);
sim_number = 9;
showFigs = 'off';
err = 2;
isScaled = true;
all_amr = false;
art_stab = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. load data (in form [level, time, q[1], x, y, aux] R:(N x 6))
disp('Loading data...')
tstart = tic;
gauges_struct = load(['matlab_matrix' num2str(sim_number) '.mat']);
gfields  = fieldnames(gauges_struct);
if sim_number >=4
    params.GC_time = gauges_struct.(gfields{1});
    gauges_struct = rmfield(gauges_struct, gfields{1});
    gfields  = fieldnames(gauges_struct);
else
    params.GC_time = 0;
end
runtimes.('LoadData') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. data analysis: adaptive mesh refinement
disp('Starting data analysis...')
tstart = tic;
amr = categorize_AMR(gauges_struct);
runtimes.('DataAnalysis') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. prep data
disp('Starting data preparation...')
tstart = tic;
%TO2D: slice (and concatenate) g_3D s.t. only use 2 observables INSIDE
%PREP_DATA
[G, scale, offset] = prep_data(gauges_struct, chosen_v, amr, isScaled);
params.cutoff_idx = round(size(G, 2)*0.8);
global cutoff
cutoff = params.cutoff_idx;
G_train = G(:, 1:params.cutoff_idx); % only use some of the data for training
runtimes.('DataPrep') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X. user variables
params.create_vid = false;
params.num_frames = 300;
params.num_pred = size(G, 2);
params.grid_res = 10;
g2_fields = fieldnames(amr.AMR2.gauge_numbers);
params.cutoff_time = gauges_struct.(g2_fields{1})(params.cutoff_idx, 2);
r = fliplr(4:4:68); %'Compute'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtime_all = [];
num_chosen = numel(chosen_v);
error_all = zeros(numel(r),num_chosen);
for kk = 1:numel(r)
    rank = r(kk);
    close all
    disp(["rank = " num2str(rank)])
    % 3. LDMD Offline: Lagrangian DMD
    disp('Starting Lagrangian DMD...')
    tstart = tic;
    [Phi, D, b, t_svd] = lDMD(G_train, rank);
    if art_stab
        [theta_D, rho_D] = cart2pol(imag(D), real(D));
        rho_D(rho_D > 1) = 1 - (rho_D(rho_D > 1) - 1);
        [im_D, re_D] = pol2cart(theta_D, rho_D);
        D = re_D + 1i*im_D;
    end
    plot_EV(D)
    save(join([pwd '\Matrices\LDMD Spectral Analysis\D' get_suffix1 '.mat'], ''), 'D')
    runtimes.('Offline') = toc(tstart) + t_svd;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. LDMD Online: start making predictions
    disp('Starting predictions...')
    tstart = tic();
    Y_pred = predict_DMD(G_train, Phi, D, b, params);
    clear b Phi D re_D im_D rho_D theta_D
    runtimes.('Online') = toc(tstart);
    rtime_all = [rtime_all; [runtimes.('Offline'), runtimes.('Online')]];
    % scale back to normal size
    global n_gauges
    if isScaled
        ii = 1;
        %TO2D: change jj to number of observables wanted
        for jj = (1:num_chosen)
            G(ii:jj*n_gauges, :) = (G(ii:jj*n_gauges, :) + offset(jj)) * scale(jj) ;
            Y_pred(ii:jj*n_gauges, :) = (Y_pred(ii:jj*n_gauges, :) + offset(jj)) * scale(jj);
            ii = jj + 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 5. plot error e = yn - yn_DMD
    disp('Starting error analysis...')
    tstart = tic;
    %TO2D: change how error_analysis returns variables (only want chosen observables)
    err_struct = error_analysis(G, Y_pred, chosen_v, amr, gauges_struct, params);
    runtimes.('ErrorAnalysis') = toc(tstart);
    error_new = zeros(1, num_chosen);
    for ii = 1:num_chosen
        error_new(ii) = mean(err_struct.(chosen_v(ii)));
    end
    error_all(kk, :) = error_new;

    %TO2D: only plot chosen observables
    %% plot rank dependent plots
    if ~isequal('Compute', r) && (rank == r(end)) && (numel(r) > 1)
        error_names = strcat(repmat("$\bar{e}_{", 1, num_chosen), chosen_v, "}$");
        % error vs. rank plot
        linestyles = {'b--o', 'g-.^', 'r:v'};
        figure('visible', 'on')
        hold on
        for ii = 1:num_chosen
            plot(r, error_all(:, ii), linestyles{ii}, 'LineWidth', 2)
        end
        hold off
%         plot(r, error_all(:, 1),'b--o', r, error_all(:, 2),'g-.^',...
%             r, error_all(:, 3), 'r:v','LineWidth', 2)
    %     grid on
        ttl = ['\textbf{Average Relative Error} $\mathbf{\bar{e}^n}$ \textbf{of FOM Wave vs. ROM Wave}'];
        if err == 2
            title({ttl
                    ['where $e^n = \frac{\Vert \xi^n - \xi_{DMD}^n \Vert_2}{\Vert \xi^n \Vert_2}$']}, 'interpreter', 'latex')
        else 
            title({ttl
                    ['where $e^n = \frac{\Vert \xi^n - \xi_{DMD}^n \Vert_2}{\max \vert \xi^n \vert}$']}, 'interpreter', 'latex')
        end
        xlabel('rank, r', 'interpreter', 'latex')
        ylabel('average error (\%)', 'interpreter', 'latex')
        legend(error_names, 'interpreter', 'latex', 'Location', 'best')
        ylim([min(error_all, [], 'all')*0.9, max(max(error_all))*1.1])
        suffix2 = get_suffix2;
        saveas(gcf, join([pwd '\Plots\Error' num2str(err) '\error_train' suffix2 '.png'], ''))
        save(join([pwd '\Matrices\Error' num2str(err) '\error_train' suffix2 '.mat'], ''), 'error_all')

        % runtime vs. rank plot
        figure('visible', 'on')
        hold on
        semilogy(r, params.GC_time*ones(size(rtime_all(:, 2))), 'k-x',...
            r, rtime_all(:, 1) + rtime_all(:, 2), 'g--o',...
            r, rtime_all(:, 1), 'b-.^',...
            r, rtime_all(:, 2), 'r:v', 'LineWidth', 2)
        hold off
        grid on
        title('Lagrangian DMD Runtime vs. Rank')
        xlabel('rank, r')
        ylabel('runtime (s)')
        legend('GeoClaw', 'LDMD Total', 'LDMD Offline', 'LDMD Online')
        suffix2 = get_suffix2;
        saveas(gcf, join([pwd '\Plots\Runtime\runtime_rank' suffix2 '.png'], ''))
        save(join([pwd '\Matrices\Runtime\rtime_all' suffix2 '.mat'], ''), 'rtime_all')

        % speed-up factor vs rank plot
        figure('visible', 'on')
        speed_up = params.GC_time./(rtime_all(:, 1) + rtime_all(:, 2));
        plot(r, speed_up, 'g--', 'LineWidth', 2)
        title({['\textbf{Speed-Up Factor}']
            ['$\frac{\textrm{GeoClaw runtime}}{\textrm{LDMD runtime}}$']}, 'interpreter', 'latex')
        xlabel('rank, r')
        ylabel('speed-up factor')
        suffix2 = get_suffix2;
        saveas(gcf, join([pwd '\Plots\Runtime\speedup_rank' suffix2 '.png'], ''))
        save(join([pwd '\Matrices\Runtime\speed_up' suffix2 '.mat'], ''), 'speed_up')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6. create video for comparison, prediction only, and actual only
    disp('Creating output videos...')
    tstart = tic;
    if params.create_vid
        params.T = min(size(Y_pred, 2), size(G, 2));
        create_comp_video('ROMvsFOM', Y_pred, G, params, amr, gauges_struct)
        create_video('ROM', Y_pred(:, 1:params.T), params, amr, gauges_struct)
        create_video('FOM', G(:, 1:params.T), params, amr, gauges_struct)
    end
    runtimes.('CreateVideos') = toc(tstart)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%7. program total run time analysis
figure  %('visible', showFigs)
rt = cell2mat(struct2cell(runtimes));
bar(rt)
title('Program Runtime Analysis')
ylabel('time (s)')
ylim([0, max(rt)*1.1])
set(gca,'xticklabel',fieldnames(runtimes));
xtickangle(45)
text([1:length(rt)], rt, num2str(rt,'%0.2f'),'HorizontalAlignment','center','VerticalAlignment','bottom')
suffix2 = get_suffix2;
saveas(gcf, join([pwd '\Plots\Runtime\program_runtime' suffix2 '.png'], ''))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
