function Y_pred = predict_DMD(g_3D, Phi, D, b, params)
    Y_pred = zeros(size(g_3D, 1), params.num_pred);
    Y_im = zeros(size(g_3D, 1), params.num_pred);
    Y_pred(:, 1) = g_3D(:, 1);
    Y_im(:, 1) = imag(g_3D(:, 1));
    for i = 2:params.num_pred
        Y = Phi * ((D .^ i) * b);
        Y_im(:, i) = imag(Y);
        Y_next = real(Y);
        Y_pred(:, i) = Y_next;
    end
    mean_Y_im = mean(Y_im(:))
    max_Y_im = max(abs(Y_im(:)))