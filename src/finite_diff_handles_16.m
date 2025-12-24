function [grad_handle, hess_handle] = finite_diff_handles_16(type_h,h)
 

    grad_handle = @(x) compute_gradient(x, h, type_h);
    hess_handle = @(x) compute_hessian(x, h, type_h);
end

function grad_approx = compute_gradient(x, step, type_h)
    n = length(x);
    grad_approx = zeros(n,1);
    h = step * ones(n,1);
    if type_h == 1 
        h = step * abs(x); 
    end  
   
    for i = 1:n-1
        grad_approx(i) = (2*i*sin(x(i))*sin(h(i)) + 4*cos(x(i))*sin(h(i))) / (2*h(i));
    end

    grad_approx(n) = (2*n*sin(x(n))*sin(h(n)) - 2*(n-1)*cos(x(n))*sin(h(n))) / (2*h(n));
end

function Hess_approx = compute_hessian(x, step, type_h)
    n = length(x);
    Hess_approx = sparse(n,n);
    
    h = step * ones(n,1);
    if type_h == 1 
        h = step * abs(x);
    end

    for i = 1:n-1
        Hess_approx(i,i) = (i*cos(x(i)) - 2*sin(x(i))) - h(i)*(i*sin(x(i)) + 2*cos(x(i)));
    end

    Hess_approx(n,n) = (n*cos(x(n)) - (n-1)*sin(x(n))) - h(n)*(n*sin(x(n)) - (n-1)*cos(x(n)));
end