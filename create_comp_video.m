function create_comp_video(ttl, Y_pred, Y, params, AMR, gauges_struct)
global showFigs
Y_pred = Y_pred(:, 1:params.T);
Y = Y(:, 1:params.T);

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
num_frames = min(min(params.num_frames, size(Y_predX_all, 2)), size(Y_X_all, 2));
[Y_predX, Y_predY, Y_predZ] = skip_times(Y_predX_all, Y_predY_all,...
    Y_predZ_all, num_frames);
[Y_X, Y_Y, Y_Z] = skip_times(Y_X_all, Y_Y_all, Y_Z_all, num_frames);
g2_fields = fieldnames(AMR.AMR2.gauge_numbers);
t = gauges_struct.(g2_fields{1})(1:AMR.AMR2.dt_final, 2);
t = t(1:params.T);
t = t(1:round(end/num_frames):end);

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
    num_grid_x = round(abs(min_Xframe - max_Xframe) * params.grid_res);
    min_Yframe = min(min(Y_predY(:, i)), min(Y_Y(:, i)));
    max_Yframe = max(max(Y_predY(:, i)), max(Y_Y(:, i)));
    num_grid_y = round(abs(min_Yframe - max_Yframe) * params.grid_res);
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
    if t(i) <= params.cutoff_time
        title({['\textbf{ROM vs. FOM Wave}']
            ['t = ' num2str(t(i))]}, 'interpreter', 'latex')
    else
        title({['\textbf{ROM vs. FOM Wave}']
            ['t = ' num2str(t(i))]}, 'Color', 'r', 'interpreter', 'latex')
    end
    ylabel('Direction of Flow (m)')
    xlabel({'Direction $\perp$ to Flow (m)'}, 'interpreter', 'latex')
    zlabel('Water Height (m)')
    legend('Interpolated Surface', 'Reduced Order Model', 'Full Order Model')
    writeVideo(v, getframe(gcf));
end
close(v)
legend('Interpolated Surface', 'Reduced Order Model', 'Full Order Model')
saveas(gcf, [pwd '\Plots\Waveform Snapshot\' ttl '-T'...
    num2str(t(i)) suffix2 '.png'])