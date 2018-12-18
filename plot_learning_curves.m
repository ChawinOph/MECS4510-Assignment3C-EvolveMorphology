%% plot the average learning curve
load 'results\learning_curve_package_3-3-3cube_p50_g400'

fittest_hist_run1 = [0 max(fit_hist_run1)];
fittest_hist_run2 = [0 max(fit_hist_run2)];

% store all runs 
all_eval_run = [n_eval_run1', n_eval_run2'];
all_fittest_hist_run = [fittest_hist_run1', fittest_hist_run2'];

[min_n_eval, indx_min_n_eval] = min(max(all_eval_run));

% resample to get evenly space
resam_n_eval = linspace(0, min_n_eval, length(all_fittest_hist_run));

EA_learning_curves = zeros(size(all_fittest_hist_run'));
for i = 1:size(EA_learning_curves, 1)
    EA_learning_curves(i, :) = interp1(all_eval_run(:, i), all_fittest_hist_run(:, i), resam_n_eval);
end

figure;
% plot all curve with the same x array
plot(resam_n_eval, EA_learning_curves); hold on;
xlim([0 max(resam_n_eval)]);

% plot the average curve
avg_EA_learning_curve = mean(EA_learning_curves);
sem_EA_learning_curve = std(EA_learning_curves)/sqrt(size(EA_learning_curves, 1));
% std_EA_learning_curve = std(EA_learning_curves);

plot(resam_n_eval, avg_EA_learning_curve);
% plotAvgWithErrorBar(resam_n_eval, avg_EA_learning_curve, std_EA_learning_curve, 10, 'r', 1.5);
[linehandle_EA, e_handle_EA] = plotAvgWithErrorBar(resam_n_eval, avg_EA_learning_curve, sem_EA_learning_curve, 10, 'b', 1.5);

title('Learning Curves');
xlabel('No. of Evaluations')
ylabel('Fitness: COM Horizontal Displacement (m)')
legend('GA', 'Random Search')
% legnd = legend([linehandle_EA, linehandle_RS, e_handle_EA, e_handle_RS],... 
%    {'Genetic Algorithm', 'Random Search',...
%     '$\pm1\sigma_{\bar{x}}$','$\pm1\sigma_{\bar{x}}$'});
legnd = legend([linehandle_EA,  e_handle_EA],... 
   {'Age-Fitness Pareto', ...
    '$\pm1\sigma_{\bar{x}}$'});
set(legnd,'Interpreter','latex')
legnd.NumColumns = 2;
grid on; grid minor