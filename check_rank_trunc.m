function [U, S, V, K_tilde, W, D, r] = check_rank_trunc(U, S, V, K_tilde, W, D, r, Y2)
global rank_trunc
if ((~(real(D(end) - D(end-1, end-1)) == 0))...
        || (~(imag(D(end) + D(end-1, end-1)) == 0)))
    r = r+1;
    rank_trunc = r;

    % 4. Compute K_tilde = U'*Y2*V*inv(S) = U'KU (K projected on dominant modes)
    K_tilde = U(:, 1:r)' * Y2 * V(:, 1:r) / S(1:r, 1:r);

    % 5. eig decomposition of K_tilde = eig decomp of K
    [W, D] = eig(K_tilde);
    D = sparse(D);
end
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);