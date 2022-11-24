function [sols, dist] = readpy(many)

% all files in the directory + save not duplicate names
files = dir('out/*_S2_*');
nSol = max(size(files));
j = 1;
for i = 1:nSol
    if i == 1
        fname(j,:) = files(i).name(1:15);
        j = j+1;
    elseif ~strcmp(files(i).name(1:15),files(i-1).name(1:15))
        fname(j,:) = files(i).name(1:15);
        j = j+1;
    end
end

if nargin == 0
    many = j-1;
end

sols = cell(many, 1);
dist = cell(many, 1);

for i = 1:many
    sols{i} = readNPY(strcat('out/', fname(i,:), '_S2_sols.npy'));
    dist{i} = readNPY(strcat('out/', fname(i,:), '_S2_dist.npy'));
end

end

%#ok<*AGROW> 