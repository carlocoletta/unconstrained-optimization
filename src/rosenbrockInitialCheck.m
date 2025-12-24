function rosenbrockInitialCheck()
    x0_a = [1.2, 1.2]';
    x0_b = [-1.2, 1]';
    [f,g,h] = rosenbrockFunction();

    generalFunctionPlot(f, [-2,2], [-1,3], 100, "Rosenbrock");

    %Nelder mead parameters
    rho = 1;
    chi = 2;
    gamma = 0.5;
    sigma = 0.5;
    kmax = 5000;
    tol = 1e-6;
    
    %NELDER MEAD
    %NB: pass x0 as a column vector, implement how to handle in
    
    disp('****Nelder mead analysis:****');
    
    t1_a = tic;
    [xbest_a, iter_a, fbest_a, roc_a] = nelderMead(f, x0_a, rho, chi, gamma, sigma, kmax, tol,0);
    time_a = toc(t1_a);

    xbest_a = xbest_a';
    norm_gradf_a = norm(g(xbest_a));
    roc_a = roc_a(end);
    
    disp('****Results for vector x0_a:****\n');
    disp('Final x:'), disp(xbest_a);
    disp(['Norma di gradf(xk): ', num2str(norm_gradf_a)]);
    disp(['f(xk): ', num2str(fbest_a)]);
    disp(['N. of Iterations x0: ', num2str(iter_a),'/',num2str(kmax)]);
    disp(['Rate of convergence: ', num2str(roc_a)]);
    disp(['Execution time x0_a: ', num2str(time_a), ' seconds']);

    
    t1_b = tic;
    [xbest_b, iter_b, fbest_b, roc_b] = nelderMead(f, x0_b, rho, chi, gamma, sigma, kmax, tol,0);
    time_b = toc(t1_b);

    xbest_b = xbest_b';
    roc_b = roc_b(end);
    norm_gradf_b = norm(gradient(xbest_b));
    
    disp('****Results for vector x0_b:****');
    disp('Final x:'), disp(xbest_b);
    disp(['Norma di gradf(xk): ', num2str(norm_gradf_b)]);
    disp(['f(xk): ', num2str(fbest_b)]);
    disp(['N. of Iterations x0: ', num2str(iter_b),'/',num2str(kmax)]);
    disp(['Rate of convergence: ', num2str(roc_b)]);

    disp(['Execution time x0_b: ', num2str(time_b), ' seconds']);


    %avg of the results
    fbest_avg = (fbest_b+fbest_a)/2;
    iter_avg = (iter_b+iter_a)/2;
    time_avg = (time_b+time_a)/2;
    norm_gradf_avg = (norm_gradf_b+norm_gradf_a)/2;
    roc_avg = (roc_a+roc_b )/ 2;

    disp(['Avg f(xk): ', num2str(fbest_avg)]);
    disp(['Avg Norma di gradf(xk): ', num2str(norm_gradf_avg)]); 
    disp(['Avg N. of Iterations x0: ', num2str(iter_avg),'/',num2str(kmax)]);
    disp(['Execution time x0_b: ', num2str(time_avg), ' seconds']);

    xbest_avg = (xbest_b+xbest_a)/2;


    T_a = table([xbest_a; xbest_b; xbest_avg], ...
          [fbest_a; fbest_b; fbest_avg], ...
          [iter_a; iter_b; iter_avg], ...
          [time_a; time_b; time_avg], ...
          [roc_a; roc_b; roc_avg], ...
          'VariableNames', {'xbest', 'fbest', 'n_iter', 'exec_time', 'roc'}, ...
          'RowNames', {'x_a^(0)', 'x_b^(0)', 'Avg'});  % Add custom row names


    % Display the table
    disp(T_a);
 

    %% MODIFIED NEWTON
    %Newton parameters
    rho = 0.5;
    c1 = 1e-4;
    btmax = 50;
    tolgrad = 1e-6; 
    tau_kmax = 100;
    kmax = 5000;
    
   
    disp('****Newton analysis:****');

    t1_a = tic;
    [xbest_a, xseq_a, fk_a, norm_gradf_a, btseq_a] = modified_Newton_beta(f, g, h, x0_a, kmax, rho, c1, btmax, tolgrad, tau_kmax);
    time_a = toc(t1_a);
    
    xbest_a = xbest_a';
    roc_a = rateOfConvergence([1 1], xseq_a);

    disp('****Results for vector x0_a:****');
    disp('Final x:'), disp(xbest_a);
    disp(['Norma di gradf(xk): ', num2str(norm_gradf_a)]);
    disp(['f(xk): ', num2str(fk_a)]);
    disp(['N. of Iterations x0_a: ', num2str(length(btseq_a)),'/',num2str(kmax)]);
    disp(['Execution time x0_b: ', num2str(time_a), ' seconds']);
    disp(['Rate of convergence: ', num2str(roc_a)]);

    
    
    
    t1_b = tic;
    [xbest_b, xseq_b, fk_b, norm_gradf_b, btseq_b] = modified_Newton_beta(f, g, h, x0_b, kmax, rho, c1, btmax, tolgrad, tau_kmax);
    time_b = toc(t1_b);
    
    xbest_b = xbest_b';
    roc_b = rateOfConvergence([1 1], xseq_b);

    disp('****Results for vector x0_b:****');
    disp('Final x:'), disp(xbest_b);
    disp(['Norma di gradf(xk): ', num2str(norm_gradf_b)]);
    disp(['f(xk): ', num2str(fk_b)]);
    disp(['N. of Iterations x0_b: ', num2str(length(btseq_b)),'/',num2str(kmax)]);
    disp(['Execution time x0_b: ', num2str(time_b), ' seconds']);
    disp(['Rate of convergence: ', num2str(roc_b)]);

        


    %avg of the results
    xbest_avg = (xbest_b+xbest_a)/2;
    fbest_avg = (fbest_b+fbest_a)/2;
    iter_avg = (iter_b+iter_a)/2;
    time_avg = (time_b+time_a)/2;
    norm_gradf_avg = (norm_gradf_b+norm_gradf_a)/2;
    roc_avg = (roc_a+roc_b )/ 2;

    disp(['Avg f(xk): ', num2str(fbest_avg)]);
    disp(['Avg Norma di gradf(xk): ', num2str(norm_gradf_avg)]); 
    disp(['Avg N. of Iterations x0: ', num2str(iter_avg),'/',num2str(kmax)]);
    disp(['Execution time x0_b: ', num2str(time_avg), ' seconds']);

    disp(['Rate of convergence: ', num2str(roc_avg)]);

    T_b = table([xbest_a; xbest_b; xbest_avg], ...
          [fbest_a; fbest_b; fbest_avg], ...
          [norm_gradf_a; norm_gradf_b; norm_gradf_avg], ...
          [iter_a; iter_b; iter_avg], ...
          [time_a; time_b; time_avg], ...
          [roc_a; roc_b; roc_avg], ...
          'VariableNames', {'xbest', 'fbest', 'grad_fk', 'n_iter', 'exec_time', 'roc'}, ...
          'RowNames', {'x_a^(0)', 'x_b^(0)', 'Avg'}); 


    % Display the table
    disp(T_b);

  
end