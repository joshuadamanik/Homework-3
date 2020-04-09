function G = hess_central_diff(X, eps, f)
    gk = grad_central_diff(X, eps, f);
    gkh_xp = grad_central_diff(X + [eps; 0], eps, f);
    gkh_yp = grad_central_diff(X + [0; eps], eps, f);
    
    gkh_xn = grad_central_diff(X - [eps; 0], eps, f);
    gkh_yn = grad_central_diff(X - [0; eps], eps, f);
    Y = 1/eps * [gkh_xp-gkh_xn, gkh_yp-gkh_yn];
    G = 0.5*[Y+Y'];
end