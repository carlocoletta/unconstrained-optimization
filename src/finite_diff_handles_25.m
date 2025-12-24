function [grad_handle, hess_handle] = finite_diff_handles_25(type_h,h)
  
    grad_handle = @(x) findiff_grad_25(x, h, type_h);
    
    hess_handle = @(x) findiff_hess_25(x, h, type_h);
end

function grad_approx = findiff_grad_25(x, h, type_h)
    n = length(x);
    grad_approx = zeros(n, 1);
    
    if mod(n, 2) == 0
        
        for k = 1:n
            switch type_h
                case 1
                    passok = h * abs(x(k));
                case 0
                    passok = h;
            end
    
            if mod(k, 2) == 0
                grad_approx(k, 1) = (-40 * passok * (10 * x(k - 1)^2 - 10 * x(k))) / (4 * passok);
            else
                grad_approx(k, 1) = (80 * x(k) * passok * (10 * x(k)^2 + 10 * passok^2 - 10 * x(k + 1)) + 4 * passok * (x(k) - 1)) / (4 * passok);
            end
        end
    else
       
        for k = 1:n
            switch type_h
                case 1
                    passok = h * abs(x(k));
                case 0
                    passok = h;
            end
    
            if k == 1
                grad_approx(1, 1) = (80 * x(k) * passok * (10 * x(k)^2 + 10 * passok^2 - 10 * x(k + 1)) + 4 * passok * (x(k) - 1) - 40 * passok * (10 * x(n)^2 - 10 * x(1))) / (4 * passok);
            elseif k == n
                grad_approx(k, 1) = (80 * x(k) * passok * (10 * x(k)^2 + 10 * passok^2 - 10 * x(1))) / (4 * passok);
            elseif mod(k, 2) == 0
                grad_approx(k, 1) = (-40 * passok * (10 * x(k - 1)^2 - 10 * x(k))) / (4 * passok);
            else
                grad_approx(k, 1) = (80 * x(k) * passok * (10 * x(k)^2 + 10 * passok^2 - 10 * x(k + 1)) + 4 * passok * (x(k) - 1)) / (4 * passok);
            end
        end
    end
end

function hessian_approx = findiff_hess_25(x, h, type_h)
    
    n = length(x);  
    
    
    i_indices = [];
    j_indices = [];
    values = [];
    cont = 1;

    
    for k = 1:n
        
        if k == 1 && mod(n, 2) == 1
            
            switch type_h
                case 0
                    hk = h;
                case 1
                    hk = h * abs(x(k));
            end

            H_kk = (40 * hk^2 * (10 * x(k)^2 - 10 * x(k + 1)) + 1400 * hk^4 + 2400 * hk^3 * x(k) + 800 * x(k)^2 * hk^2 + 2 * hk^2 + 200 * hk^2) / (2 * hk^2);
            i_indices(cont) = k;
            j_indices(cont) = k;
            values(cont) = H_kk;
            cont = cont + 1;

        elseif k == n && mod(n, 2) == 1
            
            
            switch type_h
                case 0
                    hk = h;
                case 1
                    hk = h * abs(x(k));
            end

            H_kk = (40 * hk^2 * (10 * x(k)^2 - 10 * x(1)) + 1400 * hk^4 + 2400 * hk^3 * x(k) + 800 * x(k)^2 * hk^2) / (2 * hk^2);
            i_indices(cont) = k;
            j_indices(cont) = k;
            values(cont) = H_kk;
            cont = cont + 1;

        elseif mod(k, 2) == 1
            
            
            switch type_h
                case 0
                    hk = h;
                case 1
                    hk = h * abs(x(k));
            end

            H_kk = (40 * hk^2 * (10 * x(k)^2 - 10 * x(k + 1)) + 1400 * hk^4 + 2400 * hk^3 * x(k) + 800 * x(k)^2 * hk^2 + 2 * hk^2) / (2 * hk^2);
            i_indices(cont) = k;
            j_indices(cont) = k;
            values(cont) = H_kk;
            cont = cont + 1;

        elseif mod(k, 2) == 0
            
            
            switch type_h
                case 0
                    hk = h;
                case 1
                    hk = h * abs(x(k));
            end

            H_kk = (200 * hk^2) / (2 * hk^2);
            i_indices(cont) = k;
            j_indices(cont) = k;
            values(cont) = H_kk;
            cont = cont + 1;
        end

        
        if mod(n, 2) == 1 && k == n
            
            
            switch type_h
                case 0
                    h1 = h;
                case 1
                    h1 = h * abs(x(1));
            end

            H_n1 = (20 * h1 * (-20 * x(k) * hk - 10 * hk^2)) / (2 * hk * h1);
            i_indices(cont) = n;
            j_indices(cont) = 1;
            values(cont) = H_n1;
            cont = cont + 1;

            %Symmetry
            i_indices(cont) = 1;
            j_indices(cont) = n;
            values(cont) = H_n1;
            cont = cont + 1;

        elseif mod(k, 2) == 1 && k < n
            
            switch type_h
                case 0
                    hk1 = h;
                case 1
                    hk1 = h * abs(x(k + 1));
            end

            H_k_k1 = (20 * hk1 * (-10 * hk^2 - 20 * hk * x(k))) / (2 * hk * hk1);
            i_indices(cont) = k;
            j_indices(cont) = k + 1;
            values(cont) = H_k_k1;
            cont = cont + 1;

            
            i_indices(cont) = k + 1;
            j_indices(cont) = k;
            values(cont) = H_k_k1;
            cont = cont + 1;
        end
    end

    hessian_approx = sparse(i_indices, j_indices, values, n, n);
end
