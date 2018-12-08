% Main Evolution
clc
clear
close all
%%
p = 10; % Population size
g = 1; % number of generations
s = 0.5; % selection pressure
m = 0.1; % proportion of children that get mutated
r = 0.2; % proportion of random individuals added to the population every gen
%%
genes= rand(5,9,p);
bots = MorphCube(genes);
sim = Simulator();
par_layers = zeros(s*p, s*p, g);
children = zeros(5,9,(1-s)*p);
fits = evaluate(sim, bots);

%%
tic
for i = 1:g
    % Evaluate
    
%     fits = bots.fitness;
    ages = bots.age;
    
    % Select
    [front, idx] = pareto_pick(s-0.5*r, fits, ages);
    parent_bots = bots(idx);
    parents = parents_bots.gene;
    par_fits = fits(idx);
    
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
    chidren_bots = MorphCube(children, parent_bots.ages + 1);
    
    %Add new random individuals
    rand_bots = MorphCube(rand(5,9,r*p));
    
    child_fits = evaluate(sim, children_bots);
    rand_fits = evaluate(sim, rand_bots);
    bots = [parent_bots, children_bots, rand_bots];
    fits = [par_fits, child_fits, rand_fits];
    shuffle_ind = randperm(length(bots));
    bots = bots(shuffle_ind);
    fits = fits(shuffle_ind);
end
toc
%%
disp('Done!!');
plot(par_layers(1,:,:))

bot_no = 17;
sim = Simulator(bots(bot_no));
figure;
sim.drawRobots()

% bots(bot_no).plotPDF();
bots(bot_no).plotMaterial();

[frames, K, V, COM, fitness] = sim.simulate_and_plot(bots(bot_no));

% export to video
myVideo = VideoWriter('MorphCube.avi');
myVideo.FrameRate = 25;  % Default 30
myVideo.Quality = 100;    % Default 75
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);