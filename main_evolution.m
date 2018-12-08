% Main Evolution
clc
clear
close all
%%
p = 50;
g = 10;
m = 0.02;
s = 0.5;

%%
genes= rand(5,9,p);
bots = MorphCube(genes, ones(1,p));
sim = simulator();
par_layers = zeros(s*p, s*p, g);
children = zeros(5,9,(1-s)*p);
evaluate(bots);

%%
for i = 1:g
    % Evaluate
    
    fits = bots.fitness;
    ages = bots.age;
    
    % Select
    [front, idx] = pareto_pick(s, fits, ages);
    parent_bots = bots(idx);
    parents = parents_bots.gene;
    
    %Record
    par_layers(:,:,i) = front;
    
    % Crossover
    for j = 1:2:length(parents)
        [children(i), children(i+1)] = crossoverMat(parents(i), parents(i+1));
    end
    
    %Mutation
    for j = 1:m*p
        children = mutate1(gene);
    end
    chidren_bots = MorphCube(children, parents.ages + 1);
    evaluate(children_bots);
    bots = [parent_bots, children_bots];
    shuffle_ind = randperm(length(bots));
    bots = bots(shuffle_ind);
end

%%
plot(par_layers(1,:,:))