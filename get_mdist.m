function m_dist = get_mdist(Y_predX_all, Y_predY_all)
    XY_mat = [Y_predX_all(:, end), Y_predY_all(:, end)];
    m_dist = mahal(XY_mat, XY_mat);
