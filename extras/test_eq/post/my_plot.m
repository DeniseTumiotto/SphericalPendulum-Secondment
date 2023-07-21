function my_plot(sol, data)

[n, ~] = size(sol);

for i = 1:n
    t0 = double(data{i}.t0);
    te = double(data{i}.te);
    dt = double(data{i}.dt);
    figure()
    for j = 1:50:floor((te-t0)/dt)
        R = reshape(sol{i}(1:9,j), [3,3]);
        vec_x(1:3,j) = R * [1; 0; 0];
        vec_y(1:3,j) = R * [0; 1; 0];
        vec_z(1:3,j) = R * [0; 0; 1]; 
    end
    [xS2, yS2, zS2] = sphere(360);
    surf(xS2, yS2, zS2, 'FaceAlpha', 0.1, 'EdgeColor','none'); 
    hold on
    plot3(vec_x(1, :), vec_x(2, :), vec_x(3, :), 'o', 'MarkerSize', 3)
    plot3(vec_y(1, :), vec_y(2, :), vec_y(3, :), 'o', 'MarkerSize', 3)
    plot3(vec_z(1, :), vec_z(2, :), vec_z(3, :), 'o', 'MarkerSize', 3)
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off
end
end
%#ok<*AGROW>