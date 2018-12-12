function plotFront(par_layers)
% This function plots the pareto fronts in alternating colors
% c = distinguishable_colors(size(par_layers,3));
gens = size(par_layers,3);
c = parula(gens);
figure; hold on; colormap(c)
for i = 1:gens
%     plot(par_layers(:,1,i), par_layers(:,2,i),'.', 'color', c(i,:));
    z = repmat(i, size(par_layers,1),1);
    scatter3(-par_layers(:,1,i), par_layers(:,2,i), z, 'MarkerFaceColor', c(i,:));
end
title('Fitness-Age Pareto Front through All Generations')
xlabel('Fitness'); ylabel('Age'); zlabel('Generation');
colorbar; grid on;
hold off;

% par_3d = reshape(par_layers,[],2,1); 
% z = reshape(repmat(reshape([1:gens], 1, 1, []), size(par_layers,1),1,1),[],1,1);
% scatter3(par_3d(:,1), par_3d(:,2), z);
end