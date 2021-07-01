function [g_2D, scale, offset] = prep_data(gauges_struct, chosen, AMR, isScaled)
    g1_fields = fieldnames(AMR.AMR1.gauge_numbers);
    g2_fields = fieldnames(AMR.AMR2.gauge_numbers);
    ng1 = numel(g1_fields);
    ng2 = numel(g2_fields);
    global all_amr n_gauges
    
    % build 3D matrix and extract wanted columns
    % headers are currently: [level, t, q[1], x, y, eta, aux]
    g_3D = zeros(AMR.AMR2.dt_final,...
        size(gauges_struct.(g2_fields{2}), 2),...
        n_gauges);  % time x states x gauges
    
    for ii = 1:ng2
        g_3D(:, :, ii) = gauges_struct.(g2_fields{ii})(1:AMR.AMR2.dt_final, :);
    end
    
    if all_amr
    unique_diff = unique(AMR.dt_diff);
    for ii = (ng2 + 1):(n_gauges)
        g1_mat = gauges_struct.(g1_fields{ii-ng2});
        for j = 1:numel(unique_diff)
            idx1 = find(AMR.dt_diff==unique_diff(j));
            idx3D = AMR.dt_idx(idx1);
            for k = 0:(unique_diff(j)-1) % TODO: instead off copy/paste missing data, interpolate
                g_3D(idx3D + k, :, ii) = g1_mat(idx1, :);
            end
        end
        idx3D = AMR.dt_idx(end);
        for k = 0:(AMR.AMR2.dt_final - AMR.dt_idx(end))
            g_3D(idx3D + k, :, ii) = g1_mat(end, :);
        end
    end
    end
    g_3D = g_3D(:, 4:6, :);  % extract [x, y, eta] where eta surface level (water height + bathmetry)

    % flatten along time dimension and transpose s.t. all times in single column 3D --> 2D
    g_3D = reshape(g_3D, [size(g_3D, 1), (size(g_3D, 2) * size(g_3D, 3))])';
    g_2d.x = g_3D(1:3:end, :); g_2d.y = g_3D(2:3:end, :); g_2d.h = g_3D(3:3:end, :);
    clear g_3D
    
    %TO2D: user selected variables only here! (move this to end, only rebuild g_3D at very end)
    num_chosen = numel(chosen);
    g_2D = zeros(num_chosen*n_gauges, size(g_2d.x, 2));
    ii = 1;
    for jj = (1:num_chosen)
        g_2D(ii:jj*n_gauges, :) = g_2d.(chosen(jj));
        ii = jj*n_gauges + 1;
    end
    
    % scale data
    %TO2D: user selected variables only here!
    scale = ones(1, num_chosen);
    offset = zeros(1, num_chosen);
    if isScaled
        maxx = zeros(num_chosen, 1);
        minn = zeros(num_chosen, 1);
      
        ii = 1;
        %TO2D: user selected variables only here! change jj's (update already sliced g_3D above)
        for jj = (1:num_chosen)
            maxx(jj) = max(g_2d.(chosen(jj)), [], 'all');
            minn(jj) = min(g_2d.(chosen(jj)), [], 'all');
            scale(jj) = (maxx(jj) - minn(jj)) / 2;
            offset(jj) = maxx(jj) ./ scale(jj) - 1;
            g_2D(ii:jj*n_gauges, :) = g_2D(ii:jj*n_gauges, :) / scale(jj) - offset(jj);
            ii = n_gauges*jj + 1;
        end
    end