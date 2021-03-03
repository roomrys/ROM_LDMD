% Physics Aware DMD
% note: maybe try removing data points that are stuck only when plotting?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
global sim_number showFigs
showFigs = 'off';
sim_number = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. load data (in form [level, time, q[1], x, y, aux] R:(N x 6))
disp('Loading data...')
tstart = tic;
gauges_struct = load(['matlab_matrix' num2str(sim_number) '.mat']);
gfields  = fieldnames(gauges_struct);
runtimes.('LoadData') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. data analysis: adaptive mesh refinement
disp('Starting data analysis...')
tstart = tic;
amd = mins_AMR(gauges_struct, gfields);
runtimes.('DataAnalysis') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. prep data
disp('Starting data preparation...')
tstart = tic;
g_3Da = prep_data(amd);
params.cutoff_idx = round(size(g_3Da, 2)*4/5);
global cutoff
cutoff = params.cutoff_idx;
g_3D = g_3Da(:, 1:params.cutoff_idx); % only use half the data for training
runtimes.('DataPrep') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X. user variables
params.create_vid = false;
params.num_frames = 300;
params.num_pred = size(g_3Da, 1);
params.grid_res = 10;
params.GC_time = 669.58; % load variable here instead of hard-coding
r = [44] %'Compute'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtime_all = [];
error_all = [];
for rank = r
% 3. Lagrangian DMD
disp('Starting Lagrangian DMD...')
tstart = tic;
[Phi, D, b] = lDMD(g_3D, rank);
runtimes.('LDMD') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. start making predictions
disp('Starting predictions...')
tstart = tic();
Y_pred = predict_DMD(g_3D, Phi, D, b, params);
runtimes.('Predict') = toc(tstart);
rtime_all = [rtime_all; [runtimes.('LDMD'), runtimes.('Predict')]];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. plot error e = yn - yn_DMD
disp('Starting error analysis...')
tstart = tic;
[e_tot, e_x, e_y, e_eta] = error_analysis(g_3Da, Y_pred, amd, params);
runtimes.('ErrorAnalysis') = toc(tstart);
error_all = [error_all; [max(e_tot), max(e_x), max(e_y), max(e_eta)]];

% plot rank dependent plots
if ~isequal('Compute', r) && (rank == r(end)) && (numel(r) > 1)
    % error vs. rank plot
    figure('visible', 'on')
    hold on
    plot(r, error_all(:, 1), 'k-', 'LineWidth', 2)
    plot(r, error_all(:, 2), 'b--', 'LineWidth', 2)
    plot(r, error_all(:, 3), 'g-.', 'LineWidth', 2)
    plot(r, error_all(:, 4), 'r:', 'LineWidth', 2)
    hold off
    title({['\textbf{Relative Error of Reference Wave vs. DMD Predicted Wave}']
            ['$\arg\max_n e^n = \arg\max_n \frac{\vert\vert \xi^n - \xi_{DMD}^n \vert\vert_2}{\vert\vert \xi^n \vert\vert_2}$']}, 'interpreter', 'latex')
    xlabel('rank, r', 'interpreter', 'latex')
    ylabel('max error (\%)', 'interpreter', 'latex')
    legend('e_{tot}', 'e_{\kappa}', 'e_{\gamma}', 'e_{\eta}')
    suffix2 = get_suffix2;
    saveas(gcf, [pwd '\Plots\Error\error_rank' suffix2 '.png'])

    % runtime vs. rank plot
    figure('visible', 'on')
    hold on
    plot(r, params.GC_time*ones(size(rtime_all(:, 2))), 'k-', 'LineWidth', 2)
    plot(r, rtime_all(:, 1) + rtime_all(:, 2), 'g--', 'LineWidth', 2)
    plot(r, rtime_all(:, 1), 'b-.', 'LineWidth', 2)
    plot(r, rtime_all(:, 2), 'r:', 'LineWidth', 2)
    hold off
    title('Lagrangian DMD Formulation and Prediction Runtime vs. Rank')
    xlabel('rank, r')
    ylabel('runtime (s)')
    legend('GeoClaw', 'LDMD and predict', 'LDMD', 'predict')
    suffix2 = get_suffix2;
    saveas(gcf, [pwd '\Plots\Runtime\runtime_rank' suffix2 '.png'])
    
    % speed-up factor vs rank plot
    figure('visible', 'on')
    hold on
    plot(r, params.GC_time./(rtime_all(:, 1) + rtime_all(:, 2)), 'g--', 'LineWidth', 2)
    hold off
    title({['\textbf{Speed-Up Factor}']
        ['$\frac{\textrm{GeoClaw runtime}}{\textrm{LDMD runtime}}$']}, 'interpreter', 'latex')
    xlabel('rank, r')
    ylabel('speed-up factor')
    suffix2 = get_suffix2;
    saveas(gcf, [pwd '\Plots\Runtime\speedup_rank' suffix2 '.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. create video for comparison, prediction only, and actual only
disp('Creating output videos...')
tstart = tic;
if params.create_vid
    T = min(size(Y_pred, 2), size(g_3Da, 2));
    create_comp_video('comparison', Y_pred(:, 1:T), g_3Da(:, 1:T), params.num_frames, params.grid_res)
    create_video('predict', Y_pred(:, 1:T), params.num_frames, params.grid_res)
    create_video('actual', g_3Da(:, 1:T), params.num_frames, params.grid_res)
end
runtimes.('CreateVideos') = toc(tstart)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6. program total run time analysis
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
saveas(gcf, [pwd '\Plots\Runtime\program_runtime' suffix2 '.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
function [e_tot, e_x, e_y, e_eta] = error_analysis(g_3D, Y_pred, amd, params)
    global showFigs
    % find minimum dimensions
    [row_g3D, col_g3D] = size(g_3D);
    [row_Ypred, col_Ypred] = size(Y_pred);
    min_col = min(col_g3D, col_Ypred);
    min_row = min(row_g3D, row_Ypred);
    % total error of each predicted observable
    e_tot = 100*vecnorm(g_3D(1:min_row, 1:min_col) - Y_pred(1:min_row, 1:min_col))...
        ./ vecnorm(g_3D(1:min_row, 1:min_col));
    e_x = 100*vecnorm(g_3D(1:3:min_row, 1:min_col) - Y_pred(1:3:min_row, 1:min_col))...
        ./ vecnorm(g_3D(1:min_row, 1:min_col));
    e_y = 100*vecnorm(g_3D(2:3:min_row, 1:min_col) - Y_pred(2:3:min_row, 1:min_col))...
        ./ vecnorm(g_3D(1:min_row, 1:min_col));
    e_eta = 100*vecnorm(g_3D(3:3:min_row, 1:min_col) - Y_pred(3:3:min_row, 1:min_col))...
        ./ vecnorm(g_3D(1:min_row, 1:min_col));
    g2_fields = fieldnames(amd.list_AMR2);
    t = amd.list_AMR2.(g2_fields{1})(:, 2);
    figure('visible', showFigs)
    plot(t, e_tot, '-k', 'LineWidth', 3)
    hold on
    plot(t, e_x, '--b', 'LineWidth', 2)
    plot(t, e_y, '-.g', 'LineWidth', 2)
    plot(t, e_eta, ':r', 'LineWidth', 2)
    x1 = xline(t(params.cutoff_idx), '-.');
    hold off
    ylim([0, 100])
    global rank_trunc
    title({['\textbf{Relative Error of Reference Wave vs. DMD Predicted Wave}']
        ['$e^n = \frac{\vert\vert \xi^n - \xi_{DMD}^n \vert\vert_2}{\vert\vert \xi^n \vert\vert_2} \qquad r = '...
        num2str(rank_trunc) '$']}, 'interpreter', 'latex')
    xlabel('time (s)', 'interpreter', 'latex')
    ylabel('error (\%)', 'interpreter', 'latex')
    legend('e_{tot}', 'e_{\kappa}', 'e_{\gamma}', 'e_{\eta}', 'train data cut-off')
    suffix2 = get_suffix2;
    saveas(gcf, [pwd '\Plots\Error\error' suffix2 '.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y_pred = predict_DMD(g_3D, Phi, D, b, params)
    Y_pred = zeros(size(g_3D, 1), params.num_pred);
    Y_pred(:, 1) = g_3D(:, 1);
    for i = 2:params.num_pred
        Y_next = real(Phi * ((D .^ i) * b));
        Y_pred(:, i) = Y_next;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g_3D = prep_data(amd)
    % using only amr 2 data for DMD based on data analysis results
    g2_fields = fieldnames(amd.list_AMR2);

    % build 3D matrix and extract wanted columns
    % headers are currently: [level, t, q[1], x, y, eta, aux]
    g_3D = ones(amd.mins_AMR2,...
        size(amd.list_AMR2.(g2_fields{2}), 2), numel(g2_fields));  % time x col x gauge_id
    for i = 1:numel(g2_fields)
        g_3D(:, :, i) = amd.list_AMR2.(g2_fields{i})(1:amd.mins_AMR2, :);
    end
    g_3D = g_3D(:, 4:6, :);  % extract [x, y, eta] where eta surface level (water height + bathmetry)

    % flatten along time dimension and transpose s.t. all times in single column 3D --> 2D
    g_3D = reshape(g_3D, [size(g_3D, 1), (size(g_3D, 2) * size(g_3D, 3))])';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gmin, min_size] = find_min_mat(gsize, min_size, gmin, gfields, k)
    if gsize < min_size
        gmin = gfields{k};
        min_size = gsize;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AMD = mins_AMR(gauges_struct, gfields) %AMR Mins Data
global showFigs
% Find which data contains AMR level 1, 2, or both
AMD.lastt_AMR1 = []; AMD.lastt_AMR2 = []; AMD.lastt_AMRb = [];
AMD.xy0_AMR1_all = []; AMD.xy0_AMR2_all = []; AMD.xy0_AMRb_all = [];
AMD.xyT_AMR2_all = [];
AMD.mins_AMR1 = Inf; AMD.mins_AMR2 = Inf; AMD.mins_AMRb = Inf;
AMD.gmin_AMR1 = 'None'; AMD.gmin_AMR2 = 'None'; AMD.gmin_AMRb = 'None';
% categorize AMR L1, L2, and L3 data
for i = 1:numel(gfields)
    g_mat = gauges_struct.(gfields{i});
    % multiply first column (AMR level #: 1 or 2) by ones vector
    sum_col1 = ones(1, size(g_mat, 1)) * g_mat(:, 1);
    % check if adds to num rows, double num rows, or neither
    if sum_col1 == size(g_mat, 1)
        AMD.list_AMR1.(gfields{i}) = g_mat;  % add to list
        AMD.lastt_AMR1 = [AMD.lastt_AMR1,...
            g_mat(size(g_mat, 1), 2)];
        AMD.xy0_AMR1_all = [AMD.xy0_AMR1_all, [g_mat(1, 4); g_mat(1, 5)]];  % get IC coor
        [AMD.gmin_AMR1, AMD.mins_AMR1] = ...
            find_min_mat(size(g_mat, 1), AMD.mins_AMR1,...
            AMD.gmin_AMR1, gfields, i);  % get min
    elseif sum_col1 == 2 * size(g_mat, 1)
        AMD.list_AMR2.(gfields{i}) = g_mat;
        AMD.lastt_AMR2 = [AMD.lastt_AMR2,...
            g_mat(size(g_mat, 1), 2)];
        AMD.xyT_AMR2_all = [AMD.xyT_AMR2_all, [g_mat(size(g_mat, 1), 4); g_mat(size(g_mat, 1), 5)]];
        AMD.xy0_AMR2_all = [AMD.xy0_AMR2_all, [g_mat(1, 4); g_mat(1, 5)]];
        [AMD.gmin_AMR2, AMD.mins_AMR2] =...
            find_min_mat(size(g_mat, 1), AMD.mins_AMR2,...
            AMD.gmin_AMR2, gfields, i);
    else
        AMD.list_AMRb.(gfields{i}) = g_mat;
        AMD.lastt_AMRb = [AMD.lastt_AMRb, g_mat(size(g_mat, 1), 2)];
        AMD.xy0_AMRb_all = [AMD.xy0_AMRb_all, [g_mat(1, 4); g_mat(1, 5)]];
        [AMD.gmin_AMRb, AMD.mins_AMRb] =...
            find_min_mat(size(g_mat, 1), AMD.mins_AMRb, AMD.gmin_AMRb, gfields, i);
    end
end

% % at final time step, get rid of outlier AMR2 gauges using Mahalanobis distance
% g2fields  = fieldnames(AMD.list_AMR2);
% m_dist = mahal(AMD.xyT_AMR2_all', AMD.xyT_AMR2_all');
% AMD.ptol = 1;
% mean_mdist = mean(m_dist);
% OOB_fields = {};
% for i = 1:numel(g2fields)
%     if m_dist(i) > (mean_mdist * AMD.ptol)
%         OOB_fields{end+1} = g2fields{i};
%     end
% end
% AMD.list_AMR2 = rmfield(AMD.list_AMR2, OOB_fields);
global num_gauges num_times
num_gauges = numel(fieldnames(AMD.list_AMR2));
num_times = AMD.mins_AMR2;

% plot ending times
suffix1 = get_suffix1;
figure('visible', showFigs)
histogram(AMD.lastt_AMR1, 'BinEdges', 0:0.1:10)
hold on
histogram(AMD.lastt_AMR2, 'BinEdges', 0:0.1:10)
histogram(AMD.lastt_AMRb, 'BinEdges', 0:0.1:10)
hold off
title('Lifetime of Gauges')
legend('AMR 1', 'AMR 2', 'AMR 1&2')
xlabel('time (t)')
ylabel('number gauges')
saveas(gcf, [pwd '\Plots\Gauges Lifetime\lifetime_gauges' suffix1 '.png'])
% plot data starting points with AMR1 as red and AMR2 as blue and AMRb as green
figure('visible', showFigs)
plot(AMD.xy0_AMR1_all(1, :), AMD.xy0_AMR1_all(2, :), 'r.')
hold on
plot(AMD.xy0_AMR2_all(1, :), AMD.xy0_AMR2_all(2, :), 'b*')
plot(AMD.xy0_AMRb_all(1, :), AMD.xy0_AMRb_all(2, :), 'go')
hold off
title('Different AMR Level Initial Conditions')
xlabel('x')
ylabel('y')
legend('AMR 1', 'AMR 2', 'AMR 1&2')
saveas(gcf, [pwd '\Plots\AMR IC\amr_ic' suffix1 '.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi, D, b] = lDMD(Y, r)
global showFigs
% 1. split Y matrix into Y1 and Y2
Y1 = Y(:, 1:(size(Y, 2) - 1));
Y2 = Y(:, 2:size(Y, 2));

% 2. svd of Y1
[U, S, V] = svd(Y1, 0);

if isequal('Compute', r)
% 3. rank truncated SVD, choose rank r based on SV energy
[rel_energy_Sr, sum_S] = sv_analysis(S);
for i=1:size(rel_energy_Sr, 1)
    if (S(i, i) / sum_S) < 10e-4
        r = i
        break
    end
end
end
global rank_trunc
rank_trunc = r;
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);

% 4. Compute K_tilde = U'*Y2*V*inv(S) = U'KU (K projected on dominant modes)
K_tilde = U' * Y2 * V / S;

% 5. eig decomposition of K_tilde = eig decomp of K
[W, D] = eig(K_tilde);
D = sparse(D);

% 6. eig vector K = Y2*V*inv(S)*W ~= U*W (W = eigenvector K_tilde)
Phi = U * W;  % note this is projected DMD mode not exact DMD mode

% get IC b = pinv(Phi)* x0
b = Phi \ Y(:, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rel_energy_Sr, sum_S] = sv_analysis(S)
global showFigs
suffix1 = get_suffix1;
% plot singular value decay on logarithmic scale
figure('visible', showFigs)
semilogy(diag(S) ./ S(1, 1), 'LineWidth', 3)
grid on
title('Logarithmic Singular Value Decay Relative to \sigma_1')
ylabel('$log\frac{\sigma_i}{\sigma_1}$', 'interpreter', 'latex')
xlabel('index (i)')
xlim([1 300])
saveas(gcf, [pwd '\Plots\Singular Value Decay\sv_decay' suffix1 '.png'])
% plot relative energy of singular values
sum_S = ones(1, size(S, 1)) * (diag(S));
rel_energy_Sr = cumsum((diag(S))) ./ sum_S;
figure('visible', showFigs)
plot(rel_energy_Sr,'LineWidth',2.0)
title('Relative Energy of Rank Truncated Singular Values')
ylabel({'$\sum_{i = 1:r} \sigma_i / \sum_{j = 1:m} \sigma_j$'},...
    'interpreter', 'latex')
xlabel('Rank Truncation (r)')
xlim([1, 10])
saveas(gcf, [pwd '\Plots\Singular Value Relative Energy\sv_energy'...
    suffix1 '.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_video(ttl, Y_pred, num_frames, grid_res)
global showFigs
% extract x, y, and surface height components from prediction
[Y_predX_all, Y_predY_all, Y_predZ_all] = extract_XYZ(Y_pred);

% at final time step, get rid of outlier points using Mahalanobis distance
m_dist = get_mdist(Y_predX_all, Y_predY_all);
mdist_min = 1; % mean(m_dist);
[Y_predX_all, Y_predY_all, Y_predZ_all] = extract_mdist_XYZ(Y_predX_all,...
    Y_predY_all, Y_predZ_all, m_dist, mdist_min);

% skip a few time steps for plotting/movie speed-up
num_frames = min(num_frames, size(Y_predX_all, 2));
[Y_predX, Y_predY, Y_predZ] = skip_times(Y_predX_all, Y_predY_all,...
    Y_predZ_all, num_frames);

% find min and max of Y_predZ
[min_X, max_X] = min_max(Y_predX);
[min_Z, max_Z] = min_max(Y_predZ);

% create video object
suffix2 = get_suffix2;
v = VideoWriter([pwd '\Videos\Individual\' ttl '-'...
    num2str(size(Y_predX_all, 2)) suffix2 '.avi']);
open(v);

for i = 1:min(size(Y_predX, 2), num_frames)
    % create linespace for meshgrid
    num_grid_x = round(abs(min(Y_predX(:, i)) - max(Y_predX(:, i))) * grid_res);
    num_grid_y = round(abs(min(Y_predY(:, i)) - max(Y_predY(:, i))) * grid_res);
    xlin = linspace(min(Y_predX(:, i)),max(Y_predX(:, i)),num_grid_x);
    ylin = linspace(min(Y_predY(:, i)),max(Y_predY(:, i)),num_grid_y);
    % create meshgrid
    [X, Y] = meshgrid(xlin, ylin);
    % interpolate to create surface
    f = scatteredInterpolant(Y_predX(:, i), Y_predY(:, i), Y_predZ(:, i));
    Z = f(X, Y);
    % create frames of video
    figure('visible', 'off')
    mesh(X, Y, Z)
    axis tight; hold on
    plot3(Y_predX(:, i),Y_predY(:, i),Y_predZ(:, i),'.','MarkerSize',15) %nonuniform
    hold off
    xlim([min_X max_X])
    zlim([min_Z max_Z])
    title({[ttl ' wave']})
    ylabel('Direction of Flow (m)')
    xlabel({'Direction $\perp$ to Flow (m)'}, 'interpreter', 'latex')
    zlabel('Water Height (m)')
    legend('Interpolated Surface', 'Datapoints', 'Location','best')
    writeVideo(v, getframe(gcf));
end
close(v)
title(['Final time snapshot - ' ttl])
ylabel('Direction of Flow (m)')
xlabel({'Direction $\perp$ to Flow (m)'}, 'interpreter', 'latex')
zlabel('Water Height (m)')
saveas(gcf, [pwd '\Plots\Waveform Snapshot\' ttl '-T'...
    num2str(mdist_min) suffix2 '.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_comp_video(ttl, Y_pred, Y, num_frames, grid_res)
global showFigs
% extract x, y, and surface height components from prediction
[Y_predX_all, Y_predY_all, Y_predZ_all] = extract_XYZ(Y_pred);
[Y_X_all, Y_Y_all, Y_Z_all] = extract_XYZ(Y);

% at final time step, get rid of outlier points using Mahalanobis distance
m_dist_pred = get_mdist(Y_predX_all, Y_predY_all);
m_dist_ref = get_mdist(Y_X_all, Y_Y_all);

% XY_mat = [[Y_predX_all(:, end); Y_X_all(:, end)],...
%     [Y_predY_all(:, end); Y_Y_all(:, end)]];
% m_dist = mahal(XY_mat, XY_mat);
% m_dist_pred = m_dist(1:size(Y_predX_all, 1));
% m_dist_ref = m_dist(size(Y_predX_all, 1):end);

mdist_min = 1; % mean(m_dist);
[Y_predX_all, Y_predY_all, Y_predZ_all] = extract_mdist_XYZ(Y_predX_all,...
    Y_predY_all, Y_predZ_all, m_dist_pred, mdist_min);
[Y_X_all, Y_Y_all, Y_Z_all] = extract_mdist_XYZ(Y_X_all, Y_Y_all,...
    Y_Z_all, m_dist_ref, mdist_min);

% skip a few time steps for plotting/movie speed-up
num_frames = min(min(num_frames, size(Y_predX_all, 2)), size(Y_X_all, 2));
[Y_predX, Y_predY, Y_predZ] = skip_times(Y_predX_all, Y_predY_all,...
    Y_predZ_all, num_frames);
[Y_X, Y_Y, Y_Z] = skip_times(Y_X_all, Y_Y_all, Y_Z_all, num_frames);

% find min and max of Y_predZ
[min_Xpred, max_Xpred] = min_max(Y_predX);
[min_X, max_X] = min_max(Y_X);
min_X = min(min_X, min_Xpred);
max_X = max(max_X, max_Xpred);
[min_Zpred, max_Zpred] = min_max(Y_predZ);
[min_Z, max_Z] = min_max(Y_Z);
min_Z = min(min_Z, min_Zpred);
max_Z = max(max_Z, max_Zpred);

% create video object
suffix2 = get_suffix2;
v = VideoWriter([pwd '\Videos\Comparison\' ttl '-'...
    num2str(size(Y_predX_all, 2)) suffix2 '.avi']);
open(v);

for i = 1:min(size(Y_predX, 2), num_frames)
    % create linespace for meshgrid
    min_Xframe = min(min(Y_predX(:, i)), min(Y_X(:, i)));
    max_Xframe = max(max(Y_predX(:, i)), max(Y_X(:, i)));
    num_grid_x = round(abs(min_Xframe - max_Xframe) * grid_res);
    min_Yframe = min(min(Y_predY(:, i)), min(Y_Y(:, i)));
    max_Yframe = max(max(Y_predY(:, i)), max(Y_Y(:, i)));
    num_grid_y = round(abs(min_Yframe - max_Yframe) * grid_res);
    xlin = linspace(min_Xframe,max_Xframe,num_grid_x);
    ylin = linspace(min_Yframe,max_Yframe,num_grid_y);
    % create meshgrid
    [X, Y] = meshgrid(xlin, ylin);
    % interpolate to create surface
    fpred = scatteredInterpolant(Y_predX(:, i), Y_predY(:, i), Y_predZ(:, i));
    Zpred = fpred(X, Y);
    % create frames of video
    figure('visible', 'off')
    mesh(X, Y, Zpred)
    axis tight; hold on
    plot3(Y_predX(:, i),Y_predY(:, i),Y_predZ(:, i),'.','MarkerSize',15) %nonuniform
    plot3(Y_X(:, i),Y_Y(:, i),Y_Z(:, i),'.','MarkerSize',15) %nonuniform
    hold off
    xlim([min_X max_X])
    zlim([min_Z max_Z])
    title({['Predict vs. Actual Comparison']})
    ylabel('Direction of Flow (m)')
    xlabel({'Direction $\perp$ to Flow (m)'}, 'interpreter', 'latex')
    zlabel('Water Height (m)')
    legend('Surface: Predicted', 'Data: Predicted', 'Data: Reference', 'Location','best')
    writeVideo(v, getframe(gcf));
end
close(v)
title({['Final time snapshot: ' ttl]})
ylabel('Direction of Flow (m)')
xlabel({'Direction $\perp$ to Flow (m)'}, 'interpreter', 'latex')
zlabel('Water Height (m)')
legend('Surface: Predicted', 'Data: Predicted', 'Data: Reference')
saveas(gcf, [pwd '\Plots\Waveform Snapshot\' ttl '-T'...
    num2str(mdist_min) suffix2 '.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_A, max_A] = min_max(A)
% expects matrix input, finds max and min element in matrix
min_A = min(min(A));
max_A = max(max(A));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y_x, Y_y, Y_z] = extract_XYZ(Y)
Y_x = Y(1:3:end, :);
Y_y = Y(2:3:end, :);
Y_z = Y(3:3:end, :);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y_x, Y_y, Y_z] = extract_mdist_XYZ(Y_x, Y_y, Y_z, m_dist, mdist_min)
Y_x = Y_x(m_dist < mdist_min, :);
Y_y = Y_y(m_dist < mdist_min, :);
Y_z = Y_z(m_dist < mdist_min, :);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m_dist = get_mdist(Y_predX_all, Y_predY_all)
    XY_mat = [Y_predX_all(:, end), Y_predY_all(:, end)];
    m_dist = mahal(XY_mat, XY_mat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y_x, Y_y, Y_z] = skip_times(Y_x, Y_y, Y_z, num_frames)
Y_x = Y_x(:, 1:round(end/num_frames):end);
Y_y = Y_y(:, 1:round(end/num_frames):end);
Y_z = Y_z(:, 1:round(end/num_frames):end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function suffix1 = get_suffix1
global sim_number num_gauges num_times
suffix1 = ['_s' num2str(sim_number) '_g' num2str(num_gauges) '_t'...
    num2str(num_times)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function suffix2 = get_suffix2
global rank_trunc cutoff
suffix2 = [num2str(get_suffix1) '_r' num2str(rank_trunc) '_c' num2str(cutoff)];
end