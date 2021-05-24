function [g_3D, scale, offset] = prep_data(gauges_struct, AMR, isScaled)
    % using only amr 2 data for DMD based on data analysis results
    g1_fields = fieldnames(AMR.AMR1.gauge_numbers);
    g2_fields = fieldnames(AMR.AMR2.gauge_numbers);
    ng1 = numel(g1_fields);
    ng2 = numel(g2_fields);

    % build 3D matrix and extract wanted columns
    % headers are currently: [level, t, q[1], x, y, eta, aux]
    global all_amr
    if all_amr
        n_gauges = ng1 + ng2;
    else
        n_gauges = ng1;
    end
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
    g_3DX = g_3D(1:3:end, :); g_3DY = g_3D(2:3:end, :); g_3Dh = g_3D(3:3:end, :);
    g_3D = [g_3DX; g_3DY; g_3Dh];
    
    % scale data
    if isScaled
        maxx = [max(g_3DX, [], 'all'),...
            max(g_3DY, [], 'all'),...
            max(g_3Dh, [], 'all')];
        minn = [min(g_3DX, [], 'all'),...
            min(g_3DY, [], 'all'),...
            min(g_3Dh, [], 'all')];
        
        scale = (maxx - minn) / 2;
        offset = maxx ./ scale - 1;
        
        ii = 1;
        for jj = (1:3)
            g_3D(ii:jj*n_gauges, :) = g_3D(ii:jj*n_gauges, :) / scale(jj) - offset(jj);
            ii = n_gauges*jj + 1;
        end
    else
        scale = ones(1, 3);
        offset = zeros(1, 3);
    end