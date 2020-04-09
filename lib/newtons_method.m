
function p = newtons_method(X, eps, f)
    g = grad_central_diff(X, eps, f);
    G = hess_central_diff(X, eps, f);
    p = -G\g;
end
