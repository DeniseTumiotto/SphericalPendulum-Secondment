%% MAIN FUNCTION FOR THE SOLUTION OF SPHERICAL PENDULUM
function main(method, N_TIME, k, damp, z0, T)
clearvars -except method N_TIME k damp z0 T
%% NUMERICAL PARAMETERS
% just a comment
% CHOOSE A METHOD
if nargin < 1
    method =  listdlg('PromptString',{'Choose a method'}, ...
        'ListString',{'explicit Lie Euler method', 'implicit Lie Euler Method', ...
        'implicit Lie Midpoint Rule', 'MATLAB ode45', 'Close the current run'}, 'SelectionMode', 'single', ...
         'ListSize', [200 160]);
end
switch method
    case 1
        my_method = "explicit Lie Euler method";
    case 2
        my_method = "implicit Lie Euler method";
    case 3
        my_method = "implicit midpoint rule";
    case 4
        my_method = "ode45";
    case 5
        clearvars
        return
    otherwise
        error('Method Not Implemented yet');
end

% DEFINITION OF NUMERICAL PARAMETERS AND TIME INTEGRATION
% initial and final time
% number of time steps
% time vector
% time step size
% maximal iteration steps for implicit methods
% tolerance (relative and absolute) for Netwon iteration
t0 = 0;
if nargin < 6
    prompt = {'Insert the final time'};
    dlgtitle = 'T end';
    definput = {'1'};
    dims = [1 40];
    T = inputdlg(prompt,dlgtitle,dims,definput);
    T = str2double(T{1});
end
if nargin < 2
    prompt = {'Insert a (integer) number of time steps'};
    dlgtitle = 'Number of time steps';
    definput = {'1000'};
    dims = [1 40];
    N_TIME = inputdlg(prompt,dlgtitle,dims,definput);
    N_TIME = str2double(N_TIME{1});
end
time = linspace(t0, T, N_TIME); 
dt = time(2) - time(1); disp(num2str(dt) + " time step size")

max_it = 10;
atol = 1e-10;
rtol = 1e-10;

%% PHYSICAL PARAMETERS

% LENGTH, MASS AND DAMPING OF THE PENDULUM
L = 1; 
m = 1;
if nargin < 3
    prompt = {'Insert a (non-negative) k value'};
    dlgtitle = 'k value';
    definput = {'0'};
    dims = [1 40];
    k = inputdlg(prompt,dlgtitle,dims,definput);
    k = str2double(k{1});
    prompt = {'Insert a (non-negative) damping value'};
    dlgtitle = 'Damping value';
    definput = {'0'};
    dims = [1 40];
    damp = inputdlg(prompt,dlgtitle,dims,definput);
    damp = str2double(damp{1});
end

%% SAVE PARAMETERS TO FILE
% clear not useful variable
clear prompt dlgtitle definput dims
% save the time stamp as a string.
% format: 'yyyyMMddTHHmmss'
timestamp = datestr(now,30);
filename = strcat('out/', timestamp, 'prm', '.mat');
save(filename)

%% DEFINITION OF USEFUL FUNCTIONS

% ENERGY EVALUATION
Energy_kinetic = @(q, w) 0.5 * m * cross(w, q)' * cross(w, q);
Energy_potential = @(q, w) potential(q, L, m);

% RHS OF THE SYSTEM, RESIDUAL AND JACOBIAN FOR IMPLICIT METHODS
f = @(v) fManiToAlgebra(v, damp, k); 
action = @(B, input) actionSE3(B, input);
exponentialMap = @(v) expSE3(v);
myRes = @(v0, v, h) residual(v0, v, h, f, action, exponentialMap, my_method);
myJac = @(v0, v, h) jacobian(v0, v, h, f, action, exponentialMap, my_method);

%% INITIALIZATION OF THE PROBLEM
if nargin < 5
    % CHOOSE AN INIZIALIZATION
    init =  listdlg('PromptString',{'Choose an initialization'}, ...
        'ListString',{'Random (q, w)', 'Random q, w=0', ...
        'Choose (q, w)', 'Small variation in q', ...
        'Small variation in (q, w)'}, 'SelectionMode', 'single', ...
         'ListSize', [200 160]);
    switch init
        case 1
            [q0, w0, z0] = initializeSE3();
        case 2
            [q0, w0, z0] = initializeZeroVel();
        case 3
            z0 = startingValue();
            q0 = z0(1:3);
            w0 = z0(4:6);
        case 4
            z0 = startingValue();
            [q0, w0, z0] = initializeSmallVariation(z0);
        case 5
            z0 = startingValue();
            [q0, w0, z0] = initializeSmallVariation(z0, true);
        otherwise
            error('Not valide choice');
    end
    clear init
else
    q0 = z0(1:3);
    w0 = z0(3+(1:3));
end

disp(['The initial configuration of this run is: ', newline, ...
    '[', num2str(q0(1)), ' ', num2str(q0(2)), ' ', num2str(q0(3)), ']', ' position', newline, ...
    '[', num2str(w0(1)), ' ', num2str(w0(2)), ' ', num2str(w0(3)), ']', ' angular velocity'])

qSol = zeros(3, N_TIME);
wSol = zeros(3, N_TIME);
zSol = zeros(6, N_TIME);
kinetic_energy = zeros(N_TIME, 1);
potential_energy = zeros(N_TIME, 1);

qSol(:, 1) = q0;
wSol(:, 1) = w0;
zSol(:, 1) = z0;
kinetic_energy(1) = Energy_kinetic(qSol(:, 1), wSol(:, 1));
potential_energy(1) = Energy_potential(qSol(:, 1), wSol(:, 1));

disp("Energy of this initial condition: " + ...
    num2str(kinetic_energy(1) + potential_energy(1)));

%% SOLUTION OF THE SYSTEM
for i = 1:N_TIME-1
    if method == 1
        zSol(:, i+1) = LieEuler(f, action, exponentialMap, zSol(:, i), dt);
    elseif method < 4
        zSol(:, i+1) = NewtonIt(myRes, myJac, zSol(:, i), dt, max_it, atol, rtol);
    else
        break
    end

    qSol(:, i+1) = zSol(1:3, i+1);
    wSol(:, i+1) = zSol(3+(1:3), i+1);
    kinetic_energy(i+1) = Energy_kinetic(qSol(:, i+1), wSol(:, i+1));
    potential_energy(i+1) = Energy_potential(qSol(:, i+1), wSol(:, i+1));
end

if method == 4
    tic
    options = odeset('AbsTol', atol, 'RelTol', rtol, 'Events', @myEventsFcn);
    zSol = ode45(@(t,y) S2rhs(y, damp, k), [t0, T], zSol(:, 1), options);
%     clear options
    clear options
    N_TIME = size(zSol.x, 2);
    dt = (T - t0) / N_TIME;
    time = linspace(t0, T, N_TIME);
    zSol = zSol.y;
    qSol = zSol(1:3, :);
    wSol = zSol(3+(1:3), :);
    delete(filename)
    save(filename)
%     qSol = deval(zSol, time, 1:3);
%     wSol = deval(zSol, time, 3+(1:3));
%     zSol = [qSol; wSol];
    toc
end

%% SAVE SOLUTION ON TXT FILE
% save the time stamp as a string.
% format: 'yyyyMMddTHHmmss'
filename = strcat('out/', timestamp, 'sol', '.mat');
save(filename, 'zSol')

% saveFile = strcat('out/', timestamp, 'sol', '.txt');
% fileID = fopen(saveFile, 'w');
% fprintf(fileID, '%.16f %.16f %.16f %.16f %.16f %.16f\n', zSol);
% fclose(fileID);

if nargin < 1
    %% TIME EVOLUTION OF THE SOLUTION
    plotsPlus = questdlg('Post processing analysis?', ...
        'Plots and Animations', 'Yes, all', ...
        'Yes, just plots', 'No, thank you!','Yes, all');
    
    switch plotsPlus
        case 'Yes, all'
            animation(qSol, N_TIME, dt);
            post_plots(qSol, wSol, kinetic_energy, potential_energy, time, my_method);
        case 'Yes, just plots'
            post_plots(qSol, wSol, kinetic_energy, potential_energy, time, my_method);
        case 'No, thank you!'
            disp('Thank you for using me! Have a nice day!')
    end
    
    %% MORE SIMULATIONS?
    run = questdlg('Do you want to perform a new simulation?', ...
        'New Simulation', 'Yes', 'No, thank you!','Yes');
    if strcmp(run, 'Yes')
        main
    else
        disp('It was a pleasure serving you!')
        clearvars
    end
end
end

%% USEFUL POST PROCESSING FUNCTIONS

function animation(q, steps, dstep)
        figure('Units','centimeters','Position',[10 10 20 15])
        % Create sphere surface
        [xS2, yS2, zS2] = sphere(360);
        h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
        h.EdgeColor = 'none';
        hold on
        for i = 1:steps
            % Plot pendulum at time step i
            plot3([0, q(1,i)], [0, q(2,i)], [0, q(3,i)], 'r*');
            xlabel("x")
            ylabel("y")
            zlabel("z")
            str = "Time evolution of the pendulums, T="+num2str(dstep*i);
            axis equal;
            title(str, 'FontSize', 16)
            pause(0.0000000000001);
        end
end

function post_plots(q, w, K, P, timeVec, method)
    figure('Name','Phase Space')
    plot(q(1,:), w(1,:), 'LineWidth', 2)
    hold on
    plot(q(2,:), w(2,:), 'LineWidth', 2)
    plot(q(3,:), w(3,:), 'LineWidth', 2)
    plot(q(1,1), w(1,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    plot(q(2,1), w(2,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    plot(q(3,1), w(3,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    legend('x-coordinate', 'y-coordinate', 'z-coordiante')
    xlabel('\textbf{q}', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel('$\dot{\textbf{q}}$', 'Interpreter', 'latex', 'FontSize', 16)
    title("Phase space using "+method, 'FontSize', 18)
    grid()
    
    figure('Name','Energy of the System')
    ax(1) = subplot(2, 1, 1);
    plot(timeVec, (K+P), timeVec, P, timeVec, K, 'LineWidth', 3)
    ylabel('Energy', 'FontSize', 16)
    legend('Total Energy', 'Potential Energy', 'Kinetic Energy')
    grid()
    ax(2) = subplot(2, 1, 2);
    plot(timeVec, (K+P), 'LineWidth', 3)
    ylabel('Total Energy', 'FontSize', 16)
    grid()
    linkaxes(ax, 'x');
    xlabel('Time', 'FontSize', 16)
    sgtitle("Energy of the system using "+method, 'FontSize', 18)
end

function [ value, isterminal, direction ] = myEventsFcn ( ~, y )

normq  = norm ( y(   1:3 ) );
normom = norm ( y(3+(1:3)) );

value = zeros ( 3, 1 );
value(1) = normq -  0.1;
value(2) = normq - 10.0;
value(3) = normom - 1.0e+6;

isterminal = ones ( size ( value ) );
direction = zeros ( size ( value ) );

end