function [F, grad, Hessian] = problema16()
    
    
    F = @(x) sum((1:length(x))' .* ((1 - cos(x)) + [0; sin(x(1:end-1))] - [sin(x(2:end)); 0])); 
    grad = @(x) computeGradient(x);
    Hessian = @(x) computeHessian(x);
end


function grad_val = computeGradient(x)
    n = length(x);
    grad_val=zeros(n,1);
    grad_val(1)=2*cos(x(1))+sin(x(1));
    for i=2:n-1
        grad_val(i)=i*sin(x(i))+2*cos(x(i));
    end
    grad_val(n)=n*sin(x(n))-(n-1)*cos(x(n));
end


function hess_val = computeHessian(x)
    n = length(x);

    
    main_diag = zeros(n, 1);
    main_diag(1) = -2*sin(x(1)) + cos(x(1)); 

    for i = 2:n-1
        main_diag(i) = i*cos(x(i))-2*sin(x(i)); 
    end

    main_diag(n) = n*cos(x(n)) + (n-1)*sin(x(n));

    hess_val = spdiags(main_diag, 0, n, n);
end
