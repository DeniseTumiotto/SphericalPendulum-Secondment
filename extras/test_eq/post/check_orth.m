function rslt = check_orth(sols, data)

n = size(sols);

for i = 1:max(n)
    t0 = double(data{i}.t0);
    te = double(data{i}.te);
    dt = double(data{i}.dt);
    for j = 1:floor(floor((te-t0)/dt)/10):floor((te-t0)/dt)
        rslt(i,j) = isorth(reshape(sols{i}(1:9,j),[3,3])); 
    end
end

end

%#ok<*AGROW>