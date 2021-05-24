function AMR2 = categorize_AMR(gauges_struct) %AMR Mins Data
global showFigs
gfields  = fieldnames(gauges_struct);
% Find which data contains AMR level 1, 2, or both
lastt_AMR1 = []; lastt_AMR2 = []; lastt_AMRb = [];
xy0_AMR1 = []; xy0_AMR2 = []; xy0_AMRb = [];
dt_final_AMR1 = Inf; AMR2.dt_final = Inf; dt_final_AMRb = Inf;
% categorize AMR L1, L2, and L3 data
for i = 1:numel(gfields)
    g_mat = gauges_struct.(gfields{i});
    % multiply first column (AMR level #: 1 or 2) by ones vector
    sum_col1 = ones(1, size(g_mat, 1)) * g_mat(:, 1);
    % check if adds to num rows, double num rows, or neither
    if sum_col1 == size(g_mat, 1)
        gauge_numbers_AMR1.(gfields{i}) = gfields{i};  % add to list
        lastt_AMR1 = [lastt_AMR1, g_mat(size(g_mat, 1), 2)];
        xy0_AMR1 = [xy0_AMR1, [g_mat(1, 4); g_mat(1, 5)]];  % get IC coor
        dt_final_AMR1 = find_min_mat(size(g_mat, 1), dt_final_AMR1);  % get min size
    elseif sum_col1 == (2 * size(g_mat, 1))
        AMR2.gauge_numbers.(gfields{i}) = gfields{i};
        lastt_AMR2 = [lastt_AMR2,...
            g_mat(size(g_mat, 1), 2)];
        xy0_AMR2 = [xy0_AMR2, [g_mat(1, 4); g_mat(1, 5)]];
        AMR2.dt_final = find_min_mat(size(g_mat, 1), AMR2.dt_final);
    else
        gauge_numbers_AMRb.(gfields{i}) = gfields{i};
        lastt_AMRb = [lastt_AMRb, g_mat(size(g_mat, 1), 2)];
        xy0_AMRb = [xy0_AMRb, [g_mat(1, 4); g_mat(1, 5)]];
        dt_final_AMRb = find_min_mat(size(g_mat, 1), dt_final_AMRb);
    end
end

global num_gauges num_times
num_gauges = numel(fieldnames(AMR2.gauge_numbers));
num_times = AMR2.dt_final;

% plot ending times
if numel(lastt_AMR1) > 0
    max_t = ceil(max(lastt_AMR1));
end
if (numel(lastt_AMR2) > 0)
    if lastt_AMR2 > max_t
        max_t = ceil(max(lastt_AMR2));
    end
end
if (numel(lastt_AMRb) > 0)
    if lastt_AMRb > max_t
        max_t = ceil(max(lastt_AMRb));
    end
end
suffix1 = get_suffix1;
figure('visible', showFigs)
subplot(3, 1, 1)
histogram(lastt_AMR1, 'BinEdges', 0:0.1:max_t)
title('AMR1')
ylabel('number gauges')
subplot(3, 1, 2)
histogram(lastt_AMR2, 'BinEdges', 0:0.1:max_t)
title('AMR2')
ylabel('number gauges')
subplot(3, 1, 3)
histogram(lastt_AMRb, 'BinEdges', 0:0.1:max_t)
title('AMR1&2')
xlabel('time, t (s)')
ylabel('number gauges')
sgtitle('Lifetime of Gauges')
saveas(gcf, [pwd '\Plots\Gauges Lifetime\lifetime_gauges' suffix1 '.png'])
% plot data starting points with AMR1 as red and AMR2 as blue and AMRb as green
figure('visible', showFigs)
if numel(xy0_AMR1) > 0
    plot(xy0_AMR1(1, :), xy0_AMR1(2, :), 'r.')
else
    plot(0, 10, 'r.')
end
hold on
if numel(xy0_AMR2) > 0
    plot(xy0_AMR2(1, :), xy0_AMR2(2, :), 'b*')
else
    plot(0, 10, 'b.')
end
if numel(xy0_AMRb) > 0
    plot(xy0_AMRb(1, :), xy0_AMRb(2, :), 'go')
else
    plot(0, 10, 'g.')
end
hold off
title('Different AMR Level Initial Conditions')
xlabel('x')
ylabel('y')
legend('AMR 1', 'AMR 2', 'AMR 1&2')
saveas(gcf, [pwd '\Plots\AMR IC\amr_ic' suffix1 '.png'])