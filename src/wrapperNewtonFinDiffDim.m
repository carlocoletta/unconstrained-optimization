function [fbest_avg, grad_norm_avg, iter_avg, time_avg, roc_avg, n_fails] = wrapperNewtonFinDiffDim(dim, problema, attiva_newpoints, n_newpoints, parameters, type_h, step)
    x0 = zeros(dim,1);

   %Newton parameters
    rho = parameters(1);
    c1 = parameters(2);
    btmax = parameters(3);
    tolgrad = parameters(4); 
    tau_kmax = parameters(5);
    kmax = parameters(6);


    %Variables for avgs
    fbest_avg = 0;
    grad_norm_avg = 0;
    iter_avg = 0;
    time_avg = 0;
    roc_avg = 0;
    n_fails = 0;

    
    switch problema
        case 1
            x0(1:2:end) = -1.2;
            x0(2:2:end) = 1.0; 
            f = problema1();
            [g,h] = finite_diff_handles_1(type_h, step);
        case 16
            x0(1:end) = 1;    
            f = problema16();
            [g,h] = finite_diff_handles_16(type_h, step);
        case 25
            x0(1:2:end) = -1.2;
            x0(2:2:end) = 1.0; 
            f = problema25();
            [g,h] = finite_diff_handles_25(type_h, step);
        otherwise
            fprint("Errore");
    end
    
    
    %x0
    t1 = tic;
    [xbest, fbest, grad_norm, roc, iter] = modified_Newton_beta(f, g, h, x0, kmax, rho, c1, btmax, tolgrad, tau_kmax);
    time = toc(t1);

    disp("**** x0 best ****");% disp(xbest');
    fbest_avg = fbest_avg + fbest;
    grad_norm_avg = grad_norm_avg + grad_norm;
    iter_avg = iter_avg + iter;
    if (iter>=kmax)
        n_fails = n_fails+1;
    end
    time_avg = time_avg + time;
    roc_avg = roc_avg + roc(end);


    printResults(grad_norm, fbest, iter, kmax, time, roc(end));

    if (attiva_newpoints)
        
        new_starting_points = generate_hypercube_points(x0, n_newpoints);
    
        for i=1:10
            x0 = new_starting_points(:, i);
    
            t1 = tic;
            [~, fbest, grad_norm, roc, iter] = modified_Newton_beta(f, g, h, x0, kmax, rho, c1, btmax, tolgrad, tau_kmax);
            time = toc(t1);
            
            fbest_avg = fbest_avg + fbest;
            grad_norm_avg = grad_norm_avg + grad_norm;
            iter_avg = iter_avg + iter;
            if (iter>=kmax)
                n_fails = n_fails+1;
            end
            time_avg = time_avg + time;
            %roc_avg = roc_avg + roc(end);
        
        
            printResults(grad_norm, fbest, iter, kmax, time, roc(end));
           
            
        end


        %print avg
        fbest_avg = fbest_avg/11;
        grad_norm_avg = grad_norm_avg/11;
        iter_avg = iter_avg/11;
        time_avg = time_avg/11;
        roc_avg = roc_avg/11;

        disp("*****AVG*****");
        printResults(grad_norm_avg, fbest_avg, iter_avg, kmax, time_avg, roc_avg);

    end


end