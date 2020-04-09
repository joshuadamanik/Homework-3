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

f = f3;

%% Initialization

x1 = -2;
x2 = -2;

eps_list = logspace(0,-5,10);
%eps_list = 0.1*ones(1,10);

k_list = 7*ones(1,10);
%k_list = 5:14;

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

        for n = 1:10

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