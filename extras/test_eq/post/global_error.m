function data = global_error(sol, data)

[n, ~] = size(sol);
index = 1;
old_step = double(data{1}.dt);

for i = 2:n
    if double(data{i}.dt) < old_step
        index = i;
        old_step = double(data{i}.dt);
    end
end

for i = 1:n
    if i == index
        continue
    end
    data{i}.R_err = norm(sol{i}(1:9,end)-sol{index}(1:9,end),2);
    data{i}.w_err = norm(sol{i}(10:end,end)-sol{index}(10:end,end),2);
end

end