%TO2D: if no output specified (no '=' upon calling parse), then update o
function mat = parse_xyh(matrix, chosen) % matrix col = time; top 1/3rd rows = x; 2nd 1/3rd = y; 3rd 1/3rd = h;
global n_gauges
ii = 1;
for jj = 1:numel(chosen)
    mat.(chosen(jj)) = matrix(ii:(jj*n_gauges), :);
    ii = jj*n_gauges + 1;
end
end