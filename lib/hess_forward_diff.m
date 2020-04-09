function G = hess_forward_diff(X, eps, f)
    gk = grad_central_diff(X, eps, f);
    gkh_x = grad_central_diff(X + [eps; 0], eps, f);
    gkh_y = grad_central_diff(X + [0; eps], eps, f);
    Y = 1/eps * [gkh_x-gk, gkh_y-gk];
    G = 0.5*[Y+Y'];
end