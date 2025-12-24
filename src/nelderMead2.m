function [xbest, iter, fbest, roc] = nelderMead2(f, x0, rho, chi, gamma, sigma, kmax, tol, enablePlotting)
    % The function implements the Nelder-Mead algorithm for optimizing an objective function f.
    % Inputs:
    %   - f: Objective function to minimize.
    %   - x0: Initial guess (starting point).
    %   - rho, chi, gamma, sigma: Algorithm parameters for reflection, expansion, contraction, and shrinkage.
    %   - kmax: Maximum number of iterations.
    %   - tol: Convergence tolerance (when to stop).
    %   - enablePlotting: Flag to enable or disable plotting of convergence.
    % Outputs:
    %   - xbest: Best solution found (final simplex point).
    %   - iter: Number of iterations performed.
    %   - fbest: Objective function value at the best solution.
    %   - roc: Rate of convergence (log-based).

    n = length(x0);  % The number of variables (dimensions).
    
    % Initialize the simplex matrix (n by n+1) where each column is a point.
    simplex = zeros(n, n+1);  
    simplex(:,1) = x0;  % The first point in the simplex is the starting point.

    % Generate the other n points of the simplex by perturbing the starting point.
    for i = 1:n
        ei = zeros(n,1);  % Create a unit vector.
        ei(i) = 0.05 * abs(x0(i));  % Perturb in the i-th direction.
        simplex(:, i+1) = x0 + ei;  % Add the perturbed point to the simplex.
    end

    % Sort the simplex points based on their function values.
    [fk_sorted, indices] = nelderSorter(simplex, n, f);
    x_best = simplex(:, indices(1));  % The best point is the one with the lowest function value.
    x_history = x_best;  % Track the history of the best solutions for ROC calculation.
    
    roc = [];  % Initialize Rate of Convergence (ROC) array.
    iter = 0;  % Start iteration counter.

    % Iterate until maximum iterations or convergence tolerance is met.
    while iter < kmax && (fk_sorted(n) - fk_sorted(1)) > tol
        iter = iter + 1;  % Increment iteration counter.

        % Compute the centroid (mean) of the n best points (exclude the worst).
        centroid = mean(simplex(:, indices(1:end-1)), 2);

        % Reflection step: Reflect the worst point across the centroid.
        xR = centroid + rho * (centroid - simplex(:, indices(end)));
        fxR = f(xR);  % Function value at the reflected point.

        % Check different conditions to decide the next move.
        if fxR >= fk_sorted(1) && fxR < fk_sorted(n)  
            % If the reflected point is better than the worst but not the best.
            simplex(:, indices(end)) = xR;
        elseif fxR < fk_sorted(1)  % If reflection gives the best point.
            % Expansion step: Try to expand the reflected point further.
            xE = centroid + chi * (xR - centroid);
            if f(xE) < fxR
                simplex(:, indices(end)) = xE;
            else
                simplex(:, indices(end)) = xR;  % Keep reflected point if expansion is worse.
            end
        else  % If the reflected point is worse than the worst.
            % Contraction step: Try contracting the worst point.
            if fxR > fk_sorted(n)
                xC = centroid + gamma * (simplex(:, indices(end)) - centroid);
            else
                xC = centroid + gamma * (xR - centroid);
            end
            if f(xC) < fk_sorted(end)
                simplex(:, indices(end)) = xC;
            else  % Shrink the simplex if contraction is worse.
                for i = 2:n+1
                    simplex(:, i) = simplex(:, indices(1)) + sigma * (simplex(:, i) - simplex(:, indices(1)));
                end
            end
        end

        % Re-sort the simplex points after the modification.
        [fk_sorted, indices] = nelderSorter(simplex, n, f);
        x_best = simplex(:, indices(1));  % Update best solution.
        x_history = [x_history, x_best];  % Track the history of best solutions.

        % Calculate Rate of Convergence (ROC) after at least 2 iterations.
        if iter > 2
            norm1 = norm(x_history(:, end) - x_history(:, end-1));  % Difference between current and previous best.
            norm2 = norm(x_history(:, end-1) - x_history(:, end-2));  % Difference between previous and second-to-last best.

            if norm2 > 0
                roc(end+1) = log(norm1) / log(norm2);  % Logarithmic rate of convergence.
            end
        end
        
        % Plot the convergence every 10 iterations if enabled.
        if enablePlotting && mod(iter, 10) == 0
            figure(1);
            plot(vecnorm(x_history - x_history(:,end), 2, 1), '-o', 'MarkerSize', 4, 'LineWidth', 1.5);
            xlabel('Iterations');
            ylabel('Norm Difference');
            title('Nelder-Mead Convergence (Vector-Based)');
            grid on;
            drawnow;
        end
    end

    % Output the best solution and function value.
    xbest = simplex(:, indices(1));  % Final best point.
    fbest = fk_sorted(1);  % Function value at the best point.
end
