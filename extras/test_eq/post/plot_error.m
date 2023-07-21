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
loglog(time_vec, R_err, 'LineWidth',3)
hold on
loglog(time_vec, w_err, 'LineWidth',3)
loglog(time_vec, time_vec, 'k:')

end