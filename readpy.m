function [sols, dist] = readpy(space, many, newFirst)
% Read python binary files of solutions
%
% :param space: which setting are we reading
% :param many: how many solution one wants to read
% :param newFirst: start reading from the newest result
%
% :returns: cell struct of solutions and riemannian distance

space = strtrim(space);

% all files in the directory + save not duplicate names
files = dir(strcat('out/*_', space, '*'));
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
% start reading from the newest
if nargin > 2 && newFirst && j > 2
    fname = flip(fname);
end
% set how many solution to read
if nargin < 2 || strcmp(many,'all')
    many = j-1;
end

sols = cell(many, 1);
dist = cell(many, 1);
% read
for i = 1:many
    sols{i} = readNPY(strcat('out/', fname(i,:), '_', space, '_sols.npy'));
    dist{i} = readNPY(strcat('out/', fname(i,:), '_', space, '_dist.npy'));
end

end

%#ok<*AGROW> 