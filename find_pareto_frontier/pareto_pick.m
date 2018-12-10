function [front, idx, idxB, rank] = pareto_pick(pick, fits, ages)
% This function uses fitness/age pareto layers to select the top pick%.
% Fitnesses and ages should be column vectors, pick should be between 0 and
% 1 and represents the percentage of individuals to select
front = [];
rank = zeros(1, length(fits));
idx = [];
idxB = [];
layer = 1;
while size(front,1)< pick*length(fits)
    [membership, member_value]=find_pareto_frontier([-fits, ages]);
    scatter(ages, fits); hold on
    scatter(member_value(:,2), -member_value(:,1), 'r'); hold off
    front = [front; member_value];
%     fits = fits(membership ==0);
%     ages = ages(membership ==0);
    idxB = [idxB; find(((member_value(:,1) == -fits) + (member_value(:,2) == ages))==2)];    
    fits(membership == 1) = 0;
    ages(membership ==1) = 1000;
    rank(membership ==1) = layer;
    idx = [idx; find(membership ==1)];
%     idxB = [idxB; find([-fits, ages] ==(member_value));
    
    layer = layer + 1;
end
if length(rank(rank>0)) > pick*length(fits)
    
front = front(1:pick*length(fits),:);
idx = idx(1:pick*length(fits),:);
idxB = idxB(1:pick*length(fits),:);
    
% Try to fix idxB. Maybe use intersect or ismember? Until then, just use idxA.
end