function [sols, data] = readpy(many, newFirst)
% Read python binary files of solutions
%
% :param many: how many solution one wants to read
% :param newFirst: start reading from the newest result
%
% :returns: cell struct of solutions and riemannian distance
clc

% all files in the directory + save not duplicate names
files = dir(strcat('../out/*_*'));
nSol = size(files, 1);
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

% set how many solution to read
if nargin < 1 || strcmp(many,'all')
    many = j-1;
end
% start reading from the newest
if nargin < 2
    newFirst = 0;
end
if newFirst && size(fname,1)>1
    fname = flip(fname);
end

sols = cell(many, 1);
data = cell(many, 1);

% read solutions
for i = 1:many
    sols{i} = readNPY(strcat('../out/', fname(i,:), '_sols.npy'));
    crr_data = py.open(strcat('../out/', fname(i,:), '_par.pkl'), 'rb');
    data{i} = struct(py.pickle.load(crr_data));
    crr_data.close()
end

end

%#ok<*AGROW>