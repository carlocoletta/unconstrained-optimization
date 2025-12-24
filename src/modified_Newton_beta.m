function [xk, fk, gradfk_norm, roc, k] = ...
    modified_Newton_beta(f, gradf, Hessf, x, kmax, rho, c1, btmax, tolgrad, tau_kmax)

    % Helper function for the Armijo condition in the line search.
    farmijo = @(fk, alpha, c1_gradfk_pk) fk + alpha * c1_gradfk_pk;

    % Initialization
    n = length(x);  % Length of the input (dimension of the problem)
    xk = x;          % Initialize the current point as the starting point
    fk = f(xk);      % Compute the function value at the starting point
    gradfk = gradf(xk);   % Compute the gradient at the starting point
    gradfk_norm = norm(gradfk);  % Compute the norm of the gradient
    Hessfk = Hessf(xk);       % Compute the Hessian at the starting point
    k = 0;            % Initialize the iteration counter
    
    roc = [];         % Initialize the rate of convergence (empty at first)
    x_prev = xk;      % Store the initial point for convergence tracking

    % Main optimization loop
    while k < kmax && gradfk_norm >= tolgrad
        
        % Regularization parameter to modify Hessian if necessary
        beta = 1e-3;  % Small constant for regularization
        min_diag = min(diag(Hessfk));  % Find the minimum diagonal element of the Hessian

        % If the minimum diagonal element is positive, no regularization needed
        if min_diag > 0
            tau = 0;
        else
            tau = -min_diag + beta;  % Regularize Hessian if necessary
        end
        
        % Adjust the regularization parameter tau to make Hessian positive definite
        for k_tau = 1:tau_kmax
            Bk = Hessfk + tau * speye(n);  % Add regularization term to Hessian
            [R, p] = chol(Bk);  % Cholesky decomposition of the regularized Hessian
            if p == 0, break; end  % If successful, break the loop
            tau = max(beta, 10 * tau);  % Increase tau if Cholesky fails
        end
        
        % Compute the Newton direction (using Cholesky factorization)
        y = -R' \ gradfk;  % Solve the linear system for y
        pk = R \ y;        % Solve the linear system for pk (the Newton direction)

        % Backtracking line search to find the optimal step size alpha
        alpha = 1;  % Initial step size
        c1_gradfk_pk = c1 * (gradfk' * pk);  % Compute c1 * (gradfk' * pk) for Armijo condition
        bt = 0;  % Initialize backtracking counter
        
        % Perform backtracking line search
        while bt < btmax && f(xk + alpha * pk) > farmijo(fk, alpha, c1_gradfk_pk)
            alpha = rho * alpha;  % Decrease alpha by a factor of rho
            bt = bt + 1;  % Increment backtracking counter
        end
        
        % If maximum backtracking iterations are reached and condition is not satisfied, break
        if bt == btmax && f(xk + alpha * pk) > farmijo(fk, alpha, c1_gradfk_pk)
            break;
        end

        % Update the solution with the step size alpha
        xnew = xk + alpha * pk;
        fk = f(xnew);  % Compute the function value at the new point
        gradfk = gradf(xnew);  % Compute the gradient at the new point
        gradfk_norm = norm(gradfk);  % Compute the norm of the new gradient
        Hessfk = Hessf(xnew);  % Compute the Hessian at the new point
        
        x_prev = [x_prev, xnew];  % Store the current solution for rate of convergence computation
        
        % Rate of convergence (ROC) computation
        if k > 1
            norm1 = norm(x_prev(:, end) - x_prev(:, end-1));  % Compute the difference between the last two points
            norm2 = norm(x_prev(:, end-1) - x_prev(:, end-2));  % Compute the difference between the previous two points
            if norm2 > 0
                roc(end+1) = log(norm1) / log(norm2);  % Logarithmic rate of convergence
            end
        end

        % Update the current point and iteration counter
        xk = xnew;  % Update current point
        k = k + 1;  % Increment iteration counter
    end
end
