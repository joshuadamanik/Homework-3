function p = steepest_decent(X, eps, f)
    g = grad_central_diff(X, eps, f);
    p = -g;
end