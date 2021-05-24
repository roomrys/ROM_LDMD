function [Y_x, Y_y, Y_z] = skip_times(Y_x, Y_y, Y_z, num_frames)
Y_x = Y_x(:, 1:round(end/num_frames):end);
Y_y = Y_y(:, 1:round(end/num_frames):end);
Y_z = Y_z(:, 1:round(end/num_frames):end);
