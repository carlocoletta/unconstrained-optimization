function new_points = generate_hypercube_points(x, n_vectors)
    n = length(x); 
    new_points = (2 * rand(n, n_vectors) - 1) + x(:); 
end
