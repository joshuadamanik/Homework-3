%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEWORK #3
% Joshua Julian Damanik (20194701)
% AE551 - Introduction to Optimal Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all;
addpath('lib');

%% Test functions

f1 = @(x1,x2) 0.5.*(x1-1).^2 + 10.*(x2-1).^2;
f2 = @(x1,x2) (1-x1).^2 + 100.*(x2-x1.^2).^2;
f3 = @(x1,x2) x1.^2 + 0.5.*x2.^2 + 3;

min_x(:,1) = [1, 1]';
min_x(:,2) = [1, 1]';
min_x(:,3) = [0, 0]';

f = f1;

%% Initialization

x1 = -2;
x2 = -2;

%eps_list = logspace(0,-5,10);
eps_list = 0.1*ones(1,10);

%k_list = 7*ones(1,10);
k_list = 5:14;

method_list = {'newton', 'gradient', 'rank1', 'dfp', 'bgfs'};

if isequal(f, f1)
    fname = 'f1';
elseif isequal(f, f2)
    fname = 'f2';
elseif isequal(f, f3)
    fname = 'f3';
end



for l = 1:length(method_list)
    method = method_list{l};
    
    fprintf('%s\t%s\n', fname, method);
    fprintf('eps\tk\titer\tX\tY\tError\n');
    
    for m = 1:min(length(eps_list), length(k_list))
        
        X = [x1, x2]';
        Xline = X;

        eps = eps_list(m);
        k = k_list(m);
        quasi = quasi_newton_class(eye(2));

        for n = 1:100

            if strcmp(method, 'gradient')
                %% Steepest Decent
                p = steepest_decent(X, eps, f);
            elseif strcmp(method, 'newton')
                %% Newton's method
                p = newtons_method(X, eps, f);
            elseif strcmp(method, 'rank1')
                p = quasi.rankone(X, eps, f);
            elseif strcmp(method, 'dfp')
                p = quasi.dfp(X, eps, f);
            elseif strcmp(method, 'bgfs')
                p = quasi.bgfs(X, eps, f);
            end

            %% Fibonacci search (Line search)
            [Xa, Xb] = unimodal_interval(X, eps, p, f);
            X_star = fibonacci_search(Xa, Xb, k, f);
            Xline(:,n+1) = X_star;

            %% Calculating error
            err = norm(X_star - X);
            if (err < eps * 0.01)
                break;
            end

            X = X_star;
        end

        %fprintf('Iteration #%d: (%.4f, %.4f)\n', n, X_star);


        %% Plotting data

        if isequal(f, f1)
            X_star_anal = min_x(:,1);
        elseif isequal(f, f2)
            X_star_anal = min_x(:,2);
        elseif isequal(f, f3)
            X_star_anal = min_x(:,3);
        end

        erms = norm(X_star_anal - X_star);

        fprintf('%e\t%d\t%d\t%.4f\t%.4f\t%e\n', eps, k, n, X_star, erms);
    end
end

function fib = Fib(k)
    list = [0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811];
    fib = list(k+1);
end

function X_star = fibonacci_search(Xa, Xb, k0, f)
    fib0 = Fib(k0);
    
    k = k0 - 2;
    fib_min = 0;
    fib_max = fib0;
    
    while k > 0
        del_fib = Fib(k);
        fib_range = [fib_min, fib_min + del_fib, fib_max - del_fib, fib_max];
        
        X_range = Xa + (Xb-Xa) * fib_range / fib0;
        Z_range = f(X_range(1,:), X_range(2,:));
        [~, imin] = min(Z_range);
        X_star = X_range(:,imin);
        
        if (imin < 3)
            fib_max = fib_range(3);
        else
            fib_min = fib_range(2);
        end
        k = k - 1;
    end

    X = X_star;
end

function g = grad_central_diff(X, eps, f)
    gx = (f(X(1) + eps, X(2)) - f(X(1) - eps, X(2))) / (2*eps);
    gy = (f(X(1), X(2) + eps) - f(X(1), X(2) - eps)) / (2*eps);
    g = [gx, gy]';
end

function G = hess_central_diff(X, eps, f)
    gk = grad_central_diff(X, eps, f);
    gkh_xp = grad_central_diff(X + [eps; 0], eps, f);
    gkh_yp = grad_central_diff(X + [0; eps], eps, f);
    
    gkh_xn = grad_central_diff(X - [eps; 0], eps, f);
    gkh_yn = grad_central_diff(X - [0; eps], eps, f);
    Y = 1/eps * [gkh_xp-gkh_xn, gkh_yp-gkh_yn];
    G = 0.5*[Y+Y'];
end

function p = newtons_method(X, eps, f)
    g = grad_central_diff(X, eps, f);
    G = hess_central_diff(X, eps, f);
    p = -G\g;
end

function p = steepest_decent(X, eps, f)
    g = grad_central_diff(X, eps, f);
    p = -g;
end

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