function [min_A, max_A] = min_max(A)
% expects matrix input, finds max and min element in matrix
min_A = min(min(A));
max_A = max(max(A));
