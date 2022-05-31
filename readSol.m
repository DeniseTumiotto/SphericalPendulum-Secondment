function [sol, param] = readSol(which)
% Read one solution of the spherical pendulum
% from a txt file in the folder out/
% which = 0 is the newest solution
% which = 1 is the oldest
% a negative which will take the solution starting from the newest
% a positive which will take the solution starting from the oldest

if nargin == 0
   which = 0; % saving the newest
end

% all txt file in the directory
files = dir('out/*.txt');

% save the desired file name
if which <= 0
    filename = files(end+which).name;
else
    filename = files(which).name;
end

% open the file and saving the data
fileID = fopen(strcat('out/', filename), 'r');

formatSpec = '%e';
sizeSol = [6 Inf];
sol = fscanf(fileID, formatSpec, sizeSol);

fclose(fileID);

% loading the parameters
param = load(strcat('out/', filename(1:end-7), 'prm.mat'));

end