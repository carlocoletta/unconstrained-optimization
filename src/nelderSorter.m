function  [fk_sorted, indices] = nelderSorter(simplex,n,f)

    fk = zeros(n+1,1);
    for i = 1:n+1
        fk(i) = f(simplex(:, i));
    end
    [fk_sorted, indices] = sort(fk);

end