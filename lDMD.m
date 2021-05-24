function [Phi, D, b, t_svd] = lDMD(Y, r)
global showFigs U_dmd S_dmd V_dmd
% 1. split Y matrix into Y1 and Y2
Y1 = Y(:, 1:(size(Y, 2) - 1));
Y2 = Y(:, 2:size(Y, 2));

global t_svd
if isempty(U_dmd)
    tstart = tic();
    % 2. svd of Y1
    [U_dmd, S_dmd, V_dmd] = svd(Y1, 0);

    [rel_energy_Sr, sum_S] = sv_analysis(S_dmd);
    if isequal('Compute', r)
        % 3. rank truncated SVD, choose rank r based on SV energy
        for i=1:size(rel_energy_Sr, 1)
            if (S_dmd(i, i) / sum_S) < 10e-4
                r = i
                break
            end
        end
    end
    t_svd = toc(tstart);
end
global rank_trunc
rank_trunc = r;

% 4. Compute K_tilde = U'*Y2*V*inv(S) = U'KU (K projected on dominant modes)
S_dmd_inv = diag(1 ./ diag(S_dmd(1:r, 1:r)));
K_tilde = U_dmd(:, 1:r)' * Y2 * V_dmd(:, 1:r) * S_dmd_inv;

% 5. eig decomposition of K_tilde = eig decomp of K
[W, D] = eig(K_tilde);
D = sparse(D);

% 345b. check rank truncation has complex conjugate
[U, S, V, K_tilde, W, D, r] = check_rank_trunc(U_dmd, S_dmd, V_dmd, K_tilde, W, D, r, Y2);

% 6. eig vector K = Y2*V*inv(S)*W ~= U*W (W = eigenvector K_tilde)
Phi = U * W;  % note this is projected DMD mode not exact DMD mode

% get IC b = pinv(Phi)* x0
b = Phi \ Y(:, 1);
