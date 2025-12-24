function printResults(grad_norm, fbest, iter, kmax, time, roc)
    
    disp(['f(xk): ', num2str(fbest)]);
    disp(['norm gradf(xk): ', num2str(grad_norm)]);
    disp([' N. of Iterations : ', num2str(iter),'/',num2str(kmax)]);
    disp(['Execution time: ', num2str(time), ' seconds']);
    disp(['Roc: ', num2str(roc)]);    

end