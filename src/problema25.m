function [F, gradF, HessF] = problema25()
    
    
    F = @(x) function_pb25(x);
    gradF = @(x) grad_pb25(x);
    HessF = @(x) hessian_pb25(x);
end


function val = function_pb25(x)
    n = length(x);
    val = 0;
    
    for k = 1:2:n-1
        val = val + 10 * (x(k)^2 - x(k+1))^2;
    end
    for k = 2:2:n-1
        val = val + (x(k-1) - 1)^2;
    end
    
    if mod(n,2) == 1
        val = val + (10*x(n)^2 - x(1))^2;
    else
        val = val + (x(n-1) - 1)^2;
    end
    
    val = 0.5 * val;
end


function grad = grad_pb25(x)
    n = length(x);
    grad = zeros(n,1);
    
    if mod(n,2) == 0
        grad(1:2:n-1) = 200*x(1:2:n-1).^3 - 200*x(1:2:n-1).*x(2:2:n) + x(1:2:n-1) - 1;
        grad(2:2:n) = -100*(x(1:2:n-1).^2 - x(2:2:n));
    else
        grad(1) = 200*x(1)^3 - 200*x(1)*x(2) + x(1) - 1 - 100*(x(n)^2 - x(1));
        grad(3:2:n-1) = 200*x(3:2:n-1).^3 - 200*x(3:2:n-1).*x(4:2:n) + x(3:2:n-1) - 1;
        grad(2:2:n-2) = -100*(x(1:2:n-3).^2 - x(2:2:n-2));
        grad(n) = 200*x(n).^3 - 200*x(n)*x(1) + x(n) - 1;
    end
end


function val = hessian_pb25(x)
    n = length(x);
    diags = zeros(n,5);
    
    
    if mod(n,2) == 0
        diags(2:2:n,1) = 100;
        diags(1:2:n-1,1) = 600*x(1:2:n-1).^2 - 200*x(2:2:n) + 1;
    else
        diags(1,1) = 600*x(1)^2 - 200*x(2) + 101;
        diags(2:2:n,1) = 100;
        diags(3:2:n-1,1) = 600*x(3:2:n-1).^2 - 200*x(4:2:n) + 1;
        diags(n,1) = 600*x(n).^2 - 200*x(1) + 1;
    end
    
   
    diags(1:2:n-1,3) = -200*x(1:2:n-1);
    diags(2:2:n-2,3) = 0;
    
    l
    diags(3:2:n-1,2) = 0;
    diags(2:2:n,2) = -200*x(1:2:n-1);
    
    if mod(n,2) == 1
        diags(1,5) = -200*x(n);
        diags(n,4) = -200*x(n);
    end
    
    val = spdiags(diags, [0,1,-1, n-1, -(n-1)], n, n);
end
