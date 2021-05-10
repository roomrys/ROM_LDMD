function [Y_x, Y_y, Y_z] = extract_XYZ(Y)
Y_x = Y(1:3:end, :);
Y_y = Y(2:3:end, :);
Y_z = Y(3:3:end, :);
