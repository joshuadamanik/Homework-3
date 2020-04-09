function [Xa, Xb] = unimodal_interval(X, eps, p, f)
    Xa = X;
    Xb = X;
    p = p / norm(p);

    while true
        lastX = X;
        X = X + eps * p;
        Xb = X;

        eps = 1.5*eps;

        if (f(X(1), X(2)) < f(lastX(1), lastX(2)))
            Xa = lastX;
        else
            break;
        end
    end
end