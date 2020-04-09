function g = grad_central_diff(X, eps, f)
    gx = (f(X(1) + eps, X(2)) - f(X(1) - eps, X(2))) / (2*eps);
    gy = (f(X(1), X(2) + eps) - f(X(1), X(2) - eps)) / (2*eps);
    g = [gx, gy]';
end