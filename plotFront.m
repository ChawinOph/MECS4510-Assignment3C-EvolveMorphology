function plotFront(par_layers)
% This function plots the pareto fronts in alternating colors
figure; hold on;
for i = 1:size(par_layers,3)
    scatter(par_layers(:,1,i), par_layers(:,2,i));
end
title('Fitness-Age Pareto Front through All Generations')
xlabel('fitness'); ylabel('age');
hold off;
end