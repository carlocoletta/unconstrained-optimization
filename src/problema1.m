function [F,gradF,HessF] = problema1()

    F = @(x) sum(100 * (-x(2:end) + x(1:end-1).^2).^2 + (x(1:end-1) - 1).^2);
    gradF = @(x) compute_gradient(x);
    HessF = @(x) compute_hessian(x);
    
    
    function g = compute_gradient(x)
        n = length(x);
        g = zeros(n, 1);
    
        % Gradient for the first element
        g(1) = -400 * x(1) * (x(2) - x(1)^2) + 2 * (x(1) - 1);
    
        % Gradient for the interior elements
        for i = 2:n-1
            g(i) = 200 * (x(i) - x(i-1)^2) - 400 * x(i) * (x(i+1) - x(i)^2) + 2 * (x(i) - 1);
        end
    
        % Gradient for the last element
        g(n) = 200 * (x(n) - x(n-1)^2);
    end

   function Hessf = compute_hessian(x)
        n = length(x);
        diags = zeros(n, 3);  % Matrix to store diagonal elements of Hessian
    
        % Diagonal and off-diagonal elements for the first row
        diags(1, 1) = 2 + 12*100*x(1)^2 - 4*100*x(2);
        diags(1, 2) = -4*100*x(1);
    
        % Diagonal and off-diagonal elements for the last row
        diags(n, 1) = 2*100;
        diags(n, 2) = -4*100*x(n-1);
    
        % Diagonal and off-diagonal elements for the second to last row
        diags(n-1, 3) = -4*100*x(n-1);
    
        % Loop to fill in the Hessian for interior points
        for k = 2:n-1
            diags(k, 1) = 2*100 + 12*100*x(k)^2 - 4*100*x(k+1) + 2;
            diags(k, 2) = -4*100*x(k-1);
            diags(k, 3) = -4*100*x(k);
        end
    
        % Return the sparse Hessian matrix
        Hessf = spdiags(diags, [-1, 0, 1], n, n);
    end



end
