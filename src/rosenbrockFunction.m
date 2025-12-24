function [f,g,h] = rosenbrockFunction()

    f = @(x) sum(100 * (-x(2:end) + x(1:end-1).^2).^2 + (x(1:end-1) - 1).^2);
    g = @(x) compute_gradient(x);
    h = @(x) compute_hessian(x);

    function gradf = compute_gradient(x)
        n = length(x);
        gradf = zeros(n,1);
        
        for k = 2:n-1
            gradf(k,1) = -2*100*(x(k-1)^2 - x(k)) + 2*(x(k) -1) +4*100*x(k)*(x(k)^2 - x(k+1));
        end
    
        gradf(1,1) = 2*(x(1) -1) + 4*100*x(1)*(x(1)^2 - x(2));
        gradf(n,1) = -2*100*(x(n-1)^2 - x(n)) ;
    
    end
    
    function Hessf = compute_hessian(x)
        n = length(x);
        diags = zeros(n,3);
    
        diags(1,1) = 2 + 12*100*x(1)^2 - 4*100*x(2);
        diags(n,1) = 2*100;
        diags(n-1,3) = -4*100*x(n-1);
        diags(n,2) = -4*100*x(n-1);
    
        for k = 2:n-1
           diags(k,1) = 2*100 + 12*100*x(k)^2 - 4*100*x(k+1) +2;
           diags(k-1,3) = -4*100*x(k-1);
           diags(k,2)= -4*100*x(k-1); 
        end
        
    
        Hessf = spdiags(diags, [0, +1, -1], n, n);
    
    end


end