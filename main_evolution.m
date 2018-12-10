% Main Evolution
clc
clear
close all
%%
p = 50; % Population size
g = 1; % number of generations
s = 0.5; % selection pressure
m = 0.02; % proportion of children that get mutated
r = 0.12; % proportion of random individuals added to the population every gen

%%
genes= rand(5,9,p);
bots = MorphCube(genes);
sim = Simulator();
par_layers = zeros(s*p*(1-r), 2, g);
divMat = zeros(5,9,g+1);
divMat(:,:,1) = std(genes, 0, 3);
children = zeros(5,9,(1-s)*p*(1-r));
tic
fits = evaluate(sim, bots);
toc
fit_hist = zeros(p,g+1);
fit_hist(:,1) = fits;
%%
tic
for i = 1:g
    % Evaluate
%     fits = bots.fitness;
    ages = [bots.age];
    
    % Select
    [front, idx] = pareto_pick(s-0.5*r, fits', ages');
    parent_bots = bots(idx);
    parents = [parent_bots.chromosome];
    parents = reshape(parents, 5,9,[]);
    par_fits = fits(idx);
    
    %Record
    par_layers(:,:,i) = front;
    
    % Crossover
    for j = 1:2:size(parents,3)
        [children(:,:,i), children(:,:,i+1)] = crossoverMat(parents(:,:,i), parents(:,:,i+1));
    end
    
    %Mutation
    for j = 1:m*p
        rand_idx = randi(size(children,3));
        children(rand_idx) = mutate1(genes(rand_idx));
    end
    children_bots = MorphCube(children, [parent_bots.age] + 1);
    
    %Add new random individuals
    rand_bots = MorphCube(rand(5,9,r*p));

    
    child_fits = evaluate(sim, children_bots);
    rand_fits = evaluate(sim, rand_bots);
    bots = [parent_bots, children_bots, rand_bots];
    fits = [par_fits, child_fits, rand_fits];
    shuffle_ind = randperm(length(bots));
    bots = bots(shuffle_ind);
    fits = fits(shuffle_ind);
    fit_hist(:,i+1) = fits;
    genes = cat(3, parents, children, reshape([rand_bots.chromosome],5,9,[]));
    divMat(:,:,i+1) = std(genes, 0, 3);
end
toc
%%
disp('Done!!');
plot(par_layers(1,:,:))
save('test_run3');

%%
[M,I] = max(fits);
bot_no = I;


div = sum(sum(divMat, 2),1);
figure; plot(reshape(div,1,[]));
bots(bot_no).plotPDF();
bots(bot_no).plotMaterial();

sim = Simulator(MorphCube(bots(bot_no).chromosome));
figure;
sim.drawRobots;

sim = Simulator();
[frames, K, V, COM, fitness] = sim.simulate_and_plot(MorphCube(bots(bot_no).chromosome));
tic
fitnesses = sim.evaluate(MorphCube(bots(bot_no).chromosome));
toc

% export to video
myVideo = VideoWriter('MorphCube.avi');
myVideo.FrameRate = 25;  % Default 30
myVideo.Quality = 100;    % Default 75
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);

% %% multiple robots in one video
% p_init_offset = repmat([0 0 2*0.15/2], size(sim_chrom, 3), 1) + ...
%     [1, 1, 0;       % Q1
%     -1  1  0;       % Q2
%     -1 -1  0;       % Q3
%      1 -1  0]*0.75; % Q4;
% sim = Simulator(MorphCube(sim_chrom, zeros(size(sim_chrom, 3)), 1:5, p_init_offset));
% figure;
% sim.drawRobots;
% 
% sim = Simulator();
% [frames, K, V, COM, fitness] = sim.simulate_and_plot(MorphCube(sim_chrom, zeros(size(sim_chrom, 3)), 1:5, p_init_offset));
% 
% % export to video
% myVideo = VideoWriter('MorphCubeParty.avi');
% myVideo.FrameRate = 25;  % Default 30
% myVideo.Quality = 100;    % Default 75
% open(myVideo);
% writeVideo(myVideo, frames);
% close(myVideo);
%% robot zoo
