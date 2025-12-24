function [grad_handle, hess_handle] = finite_diff_handles_1(type_h, step)

    grad_handle = @(x) compute_gradient(x, step, type_h);
    hess_handle = @(x) compute_hessian(x, step, type_h);
end

%We tried to explicitely evaluate the gradient and the hessian matrix of
%the function by hand and reporting our results to Matlab.
%Unfortunately with had issues with our code, so we ultimately decided to
%provide function handles using the general method and evaluating the
%function at each point
%We left our previous code commented at the bottom of this file



function grad_approx = compute_gradient(x, step, type_h)
    n = length(x);
    grad_approx = zeros(n, 1);
    h = step * ones(n, 1);

    if type_h == 1 
        h = step * abs(x); 
    end  

    % Gradient for the first element (finite difference with neighboring points)
    grad_approx(1) = (compute_f(x + h(1)) - compute_f(x - h(1))) / (2 * h(1));

    % Gradient for the interior elements
    for i = 2:n-1
        grad_approx(i) = (compute_f(x + h(i)) - compute_f(x - h(i))) / (2 * h(i));
    end

    % Gradient for the last element
    grad_approx(n) = (compute_f(x + h(n)) - compute_f(x - h(n))) / (2 * h(n));
end

function hess_approx = compute_hessian(x, step, type_h)
    n = length(x);
    hess_approx = sparse(n, n);
    h = step * ones(n, 1);

    if type_h == 1 
        h = step * abs(x); 
    end  

    % Diagonal elements (second derivatives)
    for i = 1:n
        hess_approx(i, i) = (compute_f(x + h(i)) - 2 * compute_f(x) + compute_f(x - h(i))) / (h(i)^2);
    end

    % Off-diagonal elements (mixed partial derivatives)
    for i = 1:n-1
        hess_approx(i, i+1) = (compute_f(x + h(i)) - compute_f(x - h(i))) / (2 * h(i));
        hess_approx(i+1, i) = hess_approx(i, i+1);  % Hessian is symmetric
    end
end

function fx = compute_f(x)
    % This function should be the objective function (Problem 1)
    % For example, the Rosenbrock function (Problem 1):
    fx = sum(100 * (x(2:end) - x(1:end-1).^2).^2 + (x(1:end-1) - 1).^2);
end


%{

function grad_approx = compute_gradient(x, step, type_h)
    n = length(x);
    grad_approx = zeros(n,1);
    h = step * ones(n,1);

    if type_h == 1 
        h = step * abs(x); 
    end  

    
    grad_approx(1) = 400*x(1)*(x(1)^2 + h(1)^2 - x(2)) + 2*(x(1)-1);

   
    for i = 2:n
        grad_approx(i) = (1 / (2 * h(i))) * ...
            (100 * (((x(i-1) + h(i-1))^2 - (x(i) + h(i)))^2 - ...
            ((x(i-1) - h(i-1))^2 - (x(i) - h(i)))^2) + ...
            (((x(i-1) + h(i-1)) - 1)^2 - ((x(i-1) - h(i-1)) - 1)^2));
    end
end



function hess_approx = compute_hessian(x, step, type_h)
    n = length(x);
    hess_approx = sparse(n,n);
    h = step * ones(n,1);

    if type_h == 1 
        h = step * abs(x); 
    end  

    %main diag
    hess_approx(1,1) = (100*(x(1) + 2*h(1))^4 - 1600*(h(1)^2*x(2)) +8*h(1)^2 -200*(x(1))^4 + 100*(x(1)-2*h(1))^2)/(h(1)^2);
    hess_approx(n,n) = 800*h(n)^2 - x(n-1)^2 + 2*x(n-1)+2*x(n)^2 - 2*x(n);
    for i=2:n-1
        hess_approx(i,i) = (800*h(i) - x(i-1)^2 + 2*x(i-1) + x(i)^2 -2*x(i) +100*(x(i)+2*h(i))^4 ...
        -1600*(h(i)^2)*x(i+1) + 8*h(i)^2 - 200*x(i)^4 +100*(x(i) - 2*h(i))^2)/(h(i)^2);
    end


    %upper diag
    hess_approx(1,2) = (-1600*h(1)*h(2)*x(1))/(4*h(1)^2);
    hess_approx(n-1, n) = (-400*x(n-1)*h(n-1))/h(n);
    
    for i=2:n-1
        hess_approx(i,i+1) = (-400*x(i)*h(i+1))/h(i);
    end

    hess_approx= max(hess_approx, hess_approx');
end
%}