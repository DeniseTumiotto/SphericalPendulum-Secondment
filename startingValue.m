function z0 = startingValue()
% This function retreives the initial values of precedent simulation
% or allows the user to insert a personalized initial condition

% all txt files in the directory
files = dir('out/*sol.mat');
[n, ~] = size(files);
initList = cell(n, 1);
listStr = cell(n+1, 1);

for i = 1:n
    filename = files(i).name;
    initList{i} = load(strcat('out/', filename));
    crrInitList = initList{i}.zSol;

    listStr{i} = strcat('[', num2str(crrInitList(1, 1)), '; ', num2str(crrInitList(2, 1)), '; ', num2str(crrInitList(3, 1)), '; ', ...
                             num2str(crrInitList(4, 1)), '; ', num2str(crrInitList(5, 1)), '; ', num2str(crrInitList(6, 1)), ']');
end

listStr{n+1} = 'Insert personal.';

init =  listdlg('PromptString',{'Choose an initialization'}, ...
    'ListString',listStr, 'SelectionMode', 'single', 'ListSize', [500 200]);

if init > n
    answer = inputdlg('Enter six space-separated numbers:', 'Initial values', [1 50]);
    z0 = str2num(answer{1});
    z0 = transpose(z0);
else
    z0 = initList{init}.zSol(:,1);
end

end