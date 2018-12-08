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
bots = MorphCube(genes);
sim = Simulator();
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