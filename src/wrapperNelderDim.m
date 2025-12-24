function [fbest_avg,iter_avg,time_avg,roc,n_failures] =wrapperNelderDim(dim, problema, f, attiva_newpoints, n_newpoints, parameters)
    x0 = zeros(dim,1);

    %Nelder mead parameters
    rho = parameters(1);
    chi = parameters(2);
    gamma = parameters(3);
    sigma = parameters(4);
    kmax = parameters(5);
    tol = parameters(6);

    %Taking x0 dependign on problem
    switch problema
        case 1
            x0(1:2:end) = -1.2;
            x0(2:2:end) = 1.0; 
        case 13
            x0(1:2:end) = -1;
            x0(2:2:end) = 1; 
        case 16
            x0(1:end) = 1; 
        case 25
            x0(1:2:end) = -1.2;
            x0(2:2:end) = 1.0;
        otherwise
            fprint("Errore");
    end
    
    %Variables for avgs
    fbest_avg = 0;
    iter_avg = 0;
    time_avg = 0;
    roc = 0;
    n_failures = 0; %I consider fail when iter == kmax


    %x0
    t1 = tic;
    [xbest, iter, fbest,roc_v] = nelderMead2(f, x0, rho, chi, gamma, sigma, kmax, tol, 0);
    time = toc(t1);
    
    %compute avgs
    fbest_avg = fbest_avg + fbest;
    iter_avg = iter_avg + iter;
    if (iter >= kmax)
        n_failures = n_failures+1;
    end
    time_avg = time_avg + time;
    roc = roc + roc_v(end);

    disp("****x0 results****");
    disp("x0 xbest:"); disp(xbest);
    printResultsNelder(fbest, iter, kmax, time, roc);

    %pertubated starts, activated through parameter
    if (attiva_newpoints)
        
        disp("****Other results****");
        %uniform distributed starting points
        new_starting_points = generate_hypercube_points(x0, n_newpoints);
    
        for i=1:10
            x0 = new_starting_points(:, i);
    
            t1 = tic;
            [~, iter, fbest, roc_v] = nelderMead2(f, x0, rho, chi, gamma, sigma, kmax, tol, 0);
            time = toc(t1);

            fbest_avg = fbest_avg + fbest;
            iter_avg = iter_avg + iter;
            if (iter >= kmax)
                 n_failures = n_failures+1;
            end
            time_avg = time_avg + time;
            roc = roc + roc_v(end);
                   
            
            printResultsNelder(fbest, iter, kmax, time, roc);

        end


        %print avg
        fbest_avg = fbest_avg/11;
        iter_avg = iter_avg/11;
        time_avg = time_avg/11;
        roc = roc/11;

        disp("****Avg results****");
        printResultsNelder(fbest, iter, kmax, time, roc);


    end


end