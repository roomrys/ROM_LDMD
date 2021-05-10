function [Y_x, Y_y, Y_z] = extract_mdist_XYZ(Y_x, Y_y, Y_z, m_dist, mdist_min)
Y_x = Y_x(m_dist < mdist_min, :);
Y_y = Y_y(m_dist < mdist_min, :);
Y_z = Y_z(m_dist < mdist_min, :);