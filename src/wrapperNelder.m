function wrapperNelder(dimensioni_nelder, problema, parameters)
    
    %Pick f dependign on problem
    switch problema
        case 1
            f = problema1();        
        case 16
            f = problema16(); 
        case 25
            f = problema25();
        otherwise
            fprint("Errore");
    end

    [fbest1, iter1, time1, roc1, fail1] = wrapperNelderDim(dimensioni_nelder(1), problema, f, 1, 10, parameters);
    [fbest2, iter2, time2, roc2, fail2] = wrapperNelderDim(dimensioni_nelder(2), problema, f, 1, 10, parameters);
    [fbest3, iter3, time3, roc3, fail3] = wrapperNelderDim(dimensioni_nelder(3), problema, f, 1, 10, parameters);

    %We print our results in a table
     T = table([fbest1; fbest2; fbest3], ...
          [iter1; iter2; iter3], ...
          [time1; time2; time3], ...
          [roc1; roc2; roc3], ...
          [fail1; fail2; fail3],...
          'VariableNames', {'fbest','n_iter', 'exec_time', 'roc', 'n_fails'}, ...
          'RowNames', {sprintf('N=%d', dimensioni_nelder(1)), ...
          sprintf('N=%d', dimensioni_nelder(2)), ...
          sprintf('N=%d', dimensioni_nelder(3))});  

    display(T);
end