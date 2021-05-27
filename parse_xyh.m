function mat = parse_xyh(matrix) % matrix col = time; top 1/3rd rows = x; 2nd 1/3rd = y; 3rd 1/3rd = h;
global n_gauges
variables = ['x', 'y', 'h'];
ii = 1;
for jj = 1:3
    mat.(variables(jj)) = matrix(ii:(jj*n_gauges), :);
    ii = jj*n_gauges + 1;
end
end