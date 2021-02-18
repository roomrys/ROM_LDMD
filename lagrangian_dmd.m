% Physics Aware DMD
% note: maybe try removing data points that are stuck only when plotting?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% user variables
create_vid = true;
num_frames = 300;
num_pred = 5000;
grid_res = 10;
global sim_number
sim_number = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. load data (in form [level, time, q[1], x, y, aux] R:(N x 6))
disp('Loading data...')
tstart = tic;
gauges_struct = load(['matlab_matrix' num2str(sim_number) '.mat']);
gfields  = fieldnames(gauges_struct);
runtimes.('LoadData') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. data analysis: adaptive mesh refinement
disp('Starting data analysis...')
tstart = tic;
amd = mins_AMR(gauges_struct, gfields);
runtimes.('DataAnalysis') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. prep data
disp('Starting data preparation...')
tstart = tic;
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
runtimes.('DataPrep') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Lagrangian DMD
disp('Starting Lagrangian DMD...')
tstart = tic;
[Phi, D] = lDMD(g_3D);
% get IC b = pinv(Phi)* x0
b = Phi \ g_3D(:, 1);
runtimes.('LDMD') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. start making predictions
disp('Starting predictions...')
tstart = tic();
Y_pred = zeros(size(g_3D, 1), num_pred);
Y_pred(:, 1) = g_3D(:, 1);
for i = 2:num_pred
    Y_next = real(Phi * diag(diag(D) .^ i) * b);
    Y_pred(:, i) = Y_next;
end
runtimes.('Predict') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. plot error e = yn - yn_DMD
disp('Starting error analysis...')
tstart = tic;
% find minimum dimensions
[row_g3D, col_g3D] = size(g_3D);
[row_Ypred, col_Ypred] = size(Y_pred);
min_col = min(col_g3D, col_Ypred);
min_row = min(row_Ypred, row_Ypred);
% total error of each predicted observable
e_tot = 100*vecnorm(g_3D(1:min_row, 1:min_col) - Y_pred(1:min_row, 1:min_col))...
    ./ vecnorm(g_3D(1:min_row, 1:min_col));
e_x = 100*vecnorm(g_3D(1:3:min_row, 1:min_col) - Y_pred(1:3:min_row, 1:min_col))...
    ./ vecnorm(g_3D(1:min_row, 1:min_col));
e_y = 100*vecnorm(g_3D(2:3:min_row, 1:min_col) - Y_pred(2:3:min_row, 1:min_col))...
    ./ vecnorm(g_3D(1:min_row, 1:min_col));
e_eta = 100*vecnorm(g_3D(3:3:min_row, 1:min_col) - Y_pred(3:3:min_row, 1:min_col))...
    ./ vecnorm(g_3D(1:min_row, 1:min_col));
figure()
plot(e_tot, '-k', 'LineWidth', 3)
hold on
plot(e_x, '--b', 'LineWidth', 2)
plot(e_y, '-.g', 'LineWidth', 2)
plot(e_eta, ':r', 'LineWidth', 2)
hold off
global rank_trunc
title({['\textbf{Relative Error of Reference Wave vs. DMD Predicted Wave}']
    ['$e^n = \Big \vert \frac{y^n - y_{DMD}^n}{y^n} \Big \vert \qquad r = '...
    num2str(rank_trunc) '$']}, 'interpreter', 'latex')
xlabel('time (s)', 'interpreter', 'latex')
ylabel('error (\%)', 'interpreter', 'latex')
legend('e_{tot}', 'e_x', 'e_y', 'e_{eta}')
xlim([0 size(e_eta, 2)])
dt = amd.list_AMR2.(g2_fields{1})(2, 2) - amd.list_AMR2.(g2_fields{1})(1, 2);
tick_loc = 0.05;
T = floor((min_col*dt)/tick_loc)*tick_loc;
xticks(0:tick_loc/dt:(T/dt))
xticklabels(string(0:tick_loc:T))
suffix2 = get_suffix2;
saveas(gcf, [pwd '\Plots\Error\error' suffix2 '.png'])
runtimes.('ErrorAnalysis') = toc(tstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. create video for comparison, prediction only, and actual only
disp('Creating output videos...')
tstart = tic;
if create_vid
    T = min(size(Y_pred, 2), size(g_3D, 2));
    create_comp_video('comparison', Y_pred(:, 1:T), g_3D(:, 1:T), num_frames, grid_res)
    create_video('predict', Y_pred(:, 1:T), num_frames, grid_res)
    create_video('actual', g_3D(:, 1:T), num_frames, grid_res)
end
runtimes.('CreateVideos') = toc(tstart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6. program total run time analysis
figure()
bar(cell2mat(struct2cell(runtimes)))
title('Program Runtime Analysis')
ylabel('time (s)')
set(gca,'xticklabel',fieldnames(runtimes));
xtickangle(45)
saveas(gcf, [pwd '\Plots\Runtime\program_runtime.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
function [gmin, min_size] = find_min_mat(gsize, min_size, gmin, gfields, k)
    if gsize < min_size
        gmin = gfields{k};
        min_size = gsize;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AMD = mins_AMR(gauges_struct, gfields) %AMR Mins Data
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
figure()
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
function [Phi, D] = lDMD(Y)
% 1. split Y matrix into Y1 and Y2
Y1 = Y(:, 1:(size(Y, 2) - 1));
Y2 = Y(:, 2:size(Y, 2));

% 2. svd of Y1
[U, S, V] = svd(Y1, 0);

% 3. rank truncated SVD, choose rank r based on SV energy
rel_energy_Sr = sv_analysis(S);
for i=1:size(rel_energy_Sr, 1)
    if rel_energy_Sr(i) >= 0.99999
        r = i
        break
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

% 6. eig vector K = Y2*V*inv(S)*W ~= U*W (W = eigenvector K_tilde)
Phi = U * W;  % note this is projected DMD mode not exact DMD mode
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rel_energy_Sr = sv_analysis(S)
suffix1 = get_suffix1;
% plot singular value decay on logarithmic scale
figure()
semilogy(diag(S) ./ S(1, 1), 'LineWidth', 3)
grid on
title('Logarithmic Singular Value Decay Relative to \sigma_1')
ylabel('$log\frac{\sigma_i}{\sigma_1}$', 'interpreter', 'latex')
xlabel('index (i)')
xlim([1 300])
saveas(gcf, [pwd '\Plots\Singular Value Decay\sv_decay' suffix1 '.png'])
% plot relative energy of singular values
sum_S = ones(1, size(S, 1)) * (diag(S) .^ 2);
rel_energy_Sr = cumsum((diag(S) .^ 2)) ./ sum_S;
figure()
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
% extract x, y, and surface height components from prediction
[Y_predX_all, Y_predY_all, Y_predZ_all] = extract_XYZ(Y_pred);

% at final time step, get rid of outlier points using Mahalanobis distance
XY_mat = [Y_predX_all(:, end), Y_predY_all(:, end)];
m_dist = mahal(XY_mat, XY_mat);
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
    writeVideo(v, getframe(gcf));
end
close(v)
title(['Final time snapshot: ' ttl])
ylabel('Direction of Flow (m)')
xlabel({'Direction $\perp$ to Flow (m)'}, 'interpreter', 'latex')
zlabel('Water Height (m)')
saveas(gcf, [pwd '\Plots\Waveform Snapshot\' ttl '-T'...
    num2str(mdist_min) suffix2 '.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_comp_video(ttl, Y_pred, Y, num_frames, grid_res)
% extract x, y, and surface height components from prediction
[Y_predX_all, Y_predY_all, Y_predZ_all] = extract_XYZ(Y_pred);
[Y_X_all, Y_Y_all, Y_Z_all] = extract_XYZ(Y);

% at final time step, get rid of outlier points using Mahalanobis distance
XY_mat = [[Y_predX_all(:, end); Y_X_all(:, end)],...
    [Y_predY_all(:, end); Y_Y_all(:, end)]];
m_dist = mahal(XY_mat, XY_mat);
m_dist_pred = m_dist(1:size(Y_predX_all, 1));
m_dist_ref = m_dist(size(Y_predX_all, 1):end);
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
suffix2 = get_suffix2
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
global rank_trunc
suffix2 = [num2str(get_suffix1) '_r' num2str(rank_trunc)];
end