function generalFunctionPlot(f, x_range, y_range, resolution, name)
    
    %Desired point range
    x1 = linspace(x_range(1), x_range(2), resolution);
    x2 = linspace(y_range(1), y_range(2), resolution);

    [X1, X2] = meshgrid(x1, x2);

    Z = arrayfun(@(a, b) f([a; b]), X1, X2);

    [~, minIndex] = min(Z(:)); 
    [row, col] = ind2sub(size(Z), minIndex); 
    x_min = X1(row, col);
    y_min = X2(row, col);
    f_min = Z(row, col);

    %Print 3d
    figure;
    surf(X1, X2, Z, 'EdgeColor', 'none');
    shading interp;
    colormap jet;
    colorbar;
    clim([min(Z(:)), min(500, max(Z(:)))]); 
    hold on;
    plot3(x_min, y_min, f_min, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

    xlabel('x1');
    ylabel('x2');
    zlabel('f(x1,x2)');
    title('3D Visualization ', name);
    view([-30, 30]);
    grid on;

    %Print 2d (top view)
    figure;
    contourf(X1, X2, Z, 50, 'LineColor', 'none'); 
    colormap jet;
    colorbar;
    clim([min(Z(:)), min(500, max(Z(:)))]);
    hold on;
    plot(x_min, y_min, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

    xlabel('x1');
    ylabel('x2');
    title('2D Visualization ', name);
    grid on;
    
end
