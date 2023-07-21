function plot_error(data)

[n, ~] = size(data);
index = 1;
old_step = double(data{1}.dt);

for i = 2:n
    if double(data{i}.dt) < old_step
        index = i;
        old_step = double(data{i}.dt);
    end
end

time_vec = zeros(n-1,1);
R_err = zeros(n-1,1);
w_err = zeros(n-1,1);

for i = 1:n
    if i == index
        continue
    end
    time_vec(i) = data{i}.dt;
    R_err(i) = data{i}.R_err;
    w_err(i) = data{i}.w_err;
end

figure()
loglog(time_vec, R_err, 'o-', 'LineWidth',3,'DisplayName','R')
hold on
loglog(time_vec, w_err, 'o-', 'LineWidth',3, 'DisplayName','\omega')
loglog(time_vec, time_vec*(R_err(end)/time_vec(end)), 'k-', 'LineWidth', 1, 'DisplayName','1^{st} order')
loglog(time_vec, time_vec.^2*(R_err(end)/time_vec(end)^2), 'k-', 'LineWidth', 1, 'DisplayName','2^{nd} order')
legend('Location','best','FontSize',16)
grid on

end