clearvars
close all
clc

[sols, params] = readUpto(2);

% % evaluate error
% error_1 = evalErr(sols(1:10), params(1:10));
% absErrq1 = error_1.q_abserr;
% absErrw1 = error_1.w_abserr;
% relErrq1 = error_1.q_relerr;
% relErrw1 = error_1.w_relerr;
% steps_1 = error_1.steps;
% 
% error_2 = evalErr(sols(11:end), params(11:end));
% absErrq2 = error_2.q_abserr;
% absErrw2 = error_2.w_abserr;
% relErrq2 = error_2.q_relerr;
% relErrw2 = error_2.w_relerr;
% steps_2 = error_2.steps;
% 
% figure()
% loglog(steps_1, absErrq1, '-or', 'LineWidth', 3, 'DisplayName','Explicit Lie Euler')
% hold on
% loglog(steps_2, absErrq2, '-ob', 'LineWidth', 3, 'DisplayName','Implicit Midpoint Rule')
% loglog(steps_1, 0.05.*steps_1, ':k', 'LineWidth', 3, 'DisplayName', 'first order')
% loglog(steps_2, 0.03.*steps_2.^2, '--k', 'LineWidth', 3, 'DisplayName', 'second order')
% legend('show', 'Location','northwest')
% xlabel('time step size', 'FontSize', 16)
% ylabel('absolute error', 'FontSize', 16)
% title('Absolute error over time step size, d=10', 'FontSize', 16)
% grid('on')


% plot distance
[m, n] = size(sols);

figure()
% for i = [m, m-5, 5, 2]
for i = 2
    distanceR = riemannianDistance(sols(i-1:i), params(i-1:i));
    distanceE = euclidDistance(sols{i-1}, sols{i});
    plot(params{i}.time, distanceR, ...
        'DisplayName', sprintf('Riemannian Distance'), ...
        'LineWidth', 1.5)
    hold on
    plot(params{i}.time, distanceE, ...
        'DisplayName', sprintf('Euclidean Distance'), ...
        'LineWidth', 1.5)
end

legend('show')
grid on