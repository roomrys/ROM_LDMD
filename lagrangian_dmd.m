% Physics Aware DMD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% user variables
create_vid = false;
num_frames = 300;
num_pred = 5000;
grid_res = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. load data (in form [level, time, q[1], x, y, aux] R:(N x 6))
gauges_struct = load('matlab_matrix2.mat');
gfields  = fieldnames(gauges_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. data analysis: adaptive mesh refinement
disp('Starting data analysis...')
tstart = tic;
amd = mins_AMR(gauges_struct, gfields);
t_data_analysis = toc(tstart)
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
t_data_prep = toc(tstart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Lagrangian DMD
disp('Starting Lagrangian DMD...')
tstart = tic;
[Phi, D] = lDMD(g_3D);
% get IC b = pinv(Phi)* x0
b = Phi \ g_3D(:, 1);
t_ldmd = toc(tstart)
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
t_predict = toc(tstart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. plot error e = yn - yn_DMD
disp('Starting error analysis...')
% find minimum dimensions
[row_g3D, col_g3D] = size(g_3D);
[row_Ypred, col_Ypred] = size(Y_pred);
min_col = min(col_g3D, col_Ypred);
min_row = min(row_Ypred, row_Ypred);
% total error of each predicted observable
e_tot = vecnorm(g_3D(1:min_row, 1:min_col) - Y_pred(1:min_row, 1:min_col));
e_x = vecnorm(g_3D(1:3:min_row, 1:min_col) - Y_pred(1:3:min_row, 1:min_col));
e_y = vecnorm(g_3D(2:3:min_row, 1:min_col) - Y_pred(2:3:min_row, 1:min_col));
e_eta = vecnorm(g_3D(3:3:min_row, 1:min_col) - Y_pred(3:3:min_row, 1:min_col));
figure()
plot(e_tot)
hold on
plot(e_x)
plot(e_y)
plot(e_eta)
hold off
title('Error: e^n = y^n - y_{DMD}^n')
xlabel('time step (n)')
ylabel('error (e)')
legend('e_{tot}', 'e_x', 'e_y', 'e_{eta}')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. create video for prediction, train, (and test)
disp('Creating output videos...')
tstart = tic;
if create_vid
    create_video('waves_predict', Y_pred, num_frames, grid_res)
    create_video('waves_train', g_3D, num_frames, grid_res)
end
t_create_videos = toc(tstart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6. program total run time analysis
t_all = [t_data_analysis, t_data_prep, t_ldmd, t_predict, t_create_videos]
figure()
bar(t_all)
title('Program Runtime Analysis')
ylabel('time (s)')
set(gca,'xticklabel',{'Data Analysis', 'Data Prep', 'lDMD', 'Predict', 'Video Creation'});
xtickangle(45)
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
AMD.lastt_AMR1 = [];
AMD.lastt_AMR2 = [];
AMD.lastt_AMRb = [];
AMD.xy0_AMR1_all = [];
AMD.xy0_AMR2_all = [];
AMD.xyT_AMR2_all = [];
AMD.xy0_AMRb_all = [];
AMD.mins_AMR1 = Inf;
AMD.mins_AMR2 = Inf;
AMD.mins_AMRb = Inf;
AMD.gmin_AMR1 = 'None';
AMD.gmin_AMR2 = 'None';
AMD.gmin_AMRb = 'None';

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

% at final time step, get rid of outlier AMR2 gauges using Mahalanobis distance
g2fields  = fieldnames(AMD.list_AMR2);
m_dist = mahal(AMD.xyT_AMR2_all', AMD.xyT_AMR2_all');
mean_mdist = mean(m_dist);
OOB_fields = {};
for i = 1:numel(g2fields)
    if m_dist(i) > mean_mdist
        OOB_fields{end+1} = g2fields{i};
    end
end
AMD.list_AMR2 = rmfield(AMD.list_AMR2, OOB_fields);

% plot ending times
histogram(AMD.lastt_AMR1, 'BinEdges', 0:0.1:10)
hold on
histogram(AMD.lastt_AMR2, 'BinEdges', 0:0.1:10)
histogram(AMD.lastt_AMRb, 'BinEdges', 0:0.1:10)
hold off
title('Lifetime of Gauges')
legend('AMR 1', 'AMR 2', 'AMR 1&2')
xlabel('time (t)')
ylabel('number gauges')

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
% plot singular value decay on logarithmic scale
figure()
semilogy(diag(S))
title('Logarithmic Singular Value Decay')
ylabel('log(\sigma_i)')
xlabel('index (i)')
% plot relative energy of singular values
sum_S = ones(1, size(S, 1)) * (diag(S) .^ 2);
rel_energy_Sr = cumsum((diag(S) .^ 2)) ./ sum_S;
figure()
plot(rel_energy_Sr,'LineWidth',2.0)
title('Relative Energy of Rank Truncated Singular Values')
ylabel({'$\sum_{i = 1:r} \sigma_i / \sum_{j = 1:m} \sigma_j$'}, 'interpreter', 'latex')
xlabel('Rank Truncation (r)')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_video(ttl, Y_pred, num_frames, grid_res)
% extract x, y, and surface height components from prediction
Y_predX_all = Y_pred(1:3:end, :);
Y_predY_all = Y_pred(2:3:end, :);
Y_predZ_all = Y_pred(3:3:end, :);

% skip a few time steps for plotting/movie speed-up
Y_predX = Y_predX_all(:, 1:round(end/num_frames):end);
Y_predY = Y_predY_all(:, 1:round(end/num_frames):end);
Y_predZ = Y_predZ_all(:, 1:round(end/num_frames):end);

% find min and max of X, Y, and Z
[min_X, max_X] = min_max(Y_predX);
[min_Y, max_Y] = min_max(Y_predY);
[min_Z, max_Z] = min_max(Y_predZ);

% create video object
v = VideoWriter([ttl '.avi']);
open(v);
% create linespace for meshgrid
num_grid_x = round(abs(min(Y_predX(:, i)) - max(Y_predX(:, i))) * grid_res);
num_grid_y = round(abs(min(Y_predY(:, i)) - max(Y_predY(:, i))) * grid_res);
%     num_grid_x = 40;
%     num_grid_y = 40;
xlin = linspace(min(Y_predX(:, i)),max(Y_predX(:, i)),num_grid_x);
ylin = linspace(min(Y_predY(:, i)),max(Y_predY(:, i)),num_grid_y);
% create meshgrid
[X, Y] = meshgrid(xlin, ylin);

for i = 1:min(size(Y_predX, 2), num_frames)
    f = scatteredInterpolant(Y_predX(:, i), Y_predY(:, i), Y_predZ(:, i));
    Z = f(X, Y);
    % create frames of video
    figure('visible', 'off')
    mesh(X, Y, Z)
    axis tight; hold on
    plot3(Y_predX(:, i),Y_predY(:, i),Y_predZ(:, i),'.','MarkerSize',15) %nonuniform
    hold off
    writeVideo(v, getframe(gcf));
end
close(v)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_A, max_A] = min_max(A)
% expects matrix input, finds max and min element in matrix
min_A = min(min(A));
max_A = max(max(A));
end