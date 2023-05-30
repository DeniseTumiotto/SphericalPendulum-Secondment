function [sols, dist, middist] = readpy(space, many, newFirst, ~, ~)
% Read python binary files of solutions
%
% :param space: which setting are we reading ('S2' or 'TS2')
% :param many: how many solution one wants to read
% :param newFirst: start reading from the newest result
%
% :returns: cell struct of solutions and riemannian distance

space = strtrim(space);

% all files in the directory + save not duplicate names
files = dir(strcat('out/*_', space, '*'));
nSol = max(size(files));
ismiddist = zeros(nSol,1);
j = 1;
for i = 1:nSol
    % check if there is 'middist'
    if strcmp(files(i).name(end-10:end),'middist.npy')
        ismiddist(i) = j-1;
    end
    if i == 1
        fname(j,:) = files(i).name(1:15);
        j = j+1;
    elseif ~strcmp(files(i).name(1:15),files(i-1).name(1:15))
        fname(j,:) = files(i).name(1:15);
        j = j+1;
    end
end
ismiddist = ismiddist(ismiddist>0);
nmid = size(ismiddist,1);
% start reading from the newest
if nargin > 2 && newFirst && j > 2
    fname = flip(fname);
    ismiddist = (j-1)*ones(nmid,1)-(ismiddist-ones(nmid,1));
end
% set how many solution to read
if nargin < 2 || strcmp(many,'all')
    many = j-1;
end

sols = cell(many, 1);
if nargout > 1
    dist = cell(many, 1);
    if nmid>0
        middist = cell(many, 1);
    else
        middist = 0;
    end
end

if nargin < 4
    % read solutions
    for i = 1:many
        sols{i} = readNPY(strcat('out/', fname(i,:), '_', space, '_sols.npy'));
        dist{i} = readNPY(strcat('out/', fname(i,:), '_', space, '_dist.npy'));
        if any(ismiddist==i)
            middist{i} = readNPY(strcat('out/', fname(i,:), '_', space, '_middist.npy'));
        end
    end
else
    % read 
    for i = 1:many
        sols{i} = readNPY(strcat('out/', fname(i,:), '_', space, '_matD.npy'));
        if nargin > 4
            dist{i} = readNPY(strcat('out/', fname(i,:), '_', space, '_matM.npy'));
        end
    end
end

end

%#ok<*AGROW> 