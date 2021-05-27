function create_video(ttl, Y_pred, params, AMR, gauges_struct)
global showFigs
Y_pred = Y_pred(:, 1:params.T);

% extract x, y, and surface height components from prediction
Y_pred_all = parse_xyh(Y_pred);
clear Y_pred
Y_predX_all = Y_pred_all.x; Y_predY_all = Y_pred_all.y; Y_predZ_all = Y_pred_all.h;
clear Y_pred_all

% at final time step, get rid of outlier points using Mahalanobis distance
m_dist = get_mdist(Y_predX_all, Y_predY_all);
mdist_min = 1; % mean(m_dist);
[Y_predX_all, Y_predY_all, Y_predZ_all] = extract_mdist_XYZ(Y_predX_all,...
    Y_predY_all, Y_predZ_all, m_dist, mdist_min);

% skip a few time steps for plotting/movie speed-up
num_frames = min(params.num_frames, size(Y_predX_all, 2));
[Y_predX, Y_predY, Y_predZ] = skip_times(Y_predX_all, Y_predY_all,...
    Y_predZ_all, num_frames);
g2_fields = fieldnames(AMR.AMR2.gauge_numbers);
t = gauges_struct.(g2_fields{1})(1:AMR.AMR2.dt_final, 2);
t = t(1:params.T);
t = t(1:round(end/num_frames):end);

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
    num_grid_x = round(abs(min(Y_predX(:, i)) - max(Y_predX(:, i))) * params.grid_res);
    num_grid_y = round(abs(min(Y_predY(:, i)) - max(Y_predY(:, i))) * params.grid_res);
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
    if t(i) <= params.cutoff_time
        title({['\textbf{' ttl ' Wave}']
            ['t = ' num2str(t(i))]}, 'interpreter', 'latex')
    else
        title({['\textbf{' ttl ' Wave}']
            ['t = ' num2str(t(i))]}, 'Color', 'r', 'interpreter', 'latex')
    end
    ylabel('Direction of Flow (m)')
    xlabel({'Direction $\perp$ to Flow (m)'}, 'interpreter', 'latex')
    zlabel('Water Height (m)')
    legend('Interpolated Surface', 'Datapoints', 'Location','best')
    writeVideo(v, getframe(gcf));
end
close(v)
saveas(gcf, [pwd '\Plots\Waveform Snapshot\' ttl '-T'...
num2str(t(i)) suffix2 '.png'])