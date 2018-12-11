% Main Evolution
clc
clear
close all
%%
p = 25; % Population size
g = 3; % number of generations

s = 0.5; % selection pressure
m = 0.02; % proportion of children that get mutated
r = 0.12; % proportion of random individuals added to the population every gen

disp(['Estimated run time: ' num2str(floor(3.5*p*(g + 1)/3600/2)) ' hr. ' num2str(round(mod(3.5*p*(g + 1)/60/2, 60))) ' mins'])

%%
n_eval = [0];
genes= rand(5,9,p);
bots = MorphCube(genes);
sim = Simulator();
par_layers = zeros(s*floor(p*(1-r)), 2, g);
divMat = zeros(5,9,g+1);
divMat(:,:,1) = std(genes, 0, 3);
children = zeros(5,9,(1-s)*floor(p*(1-r)));
tic
[fits, n_eval_gen] = sim.evaluate(bots);
n_eval = [n_eval n_eval_gen];
toc
fit_hist = zeros(p,g+1);
fit_hist(:,1) = fits;
%% start the EA
tic
for i = 1:g
    % Evaluate
%     fits = bots.fitness;
    ages = [bots.age];
    
    % Select git
    [front, idx] = pareto_pick(s-0.5*r, fits', ages');
    parent_bots = bots(idx);
%     parents = [parent_bots.chromosome];
%     parents = reshape(parents, 5,9,[]);
    par_fits = fits(idx);
    
    
    
    %Shuffle
    shuffle_ind = randperm(length(parent_bots));
    par_fits = fits(shuffle_ind);
    parent_bots = parent_bots(shuffle_ind);
    parents = [parent_bots.chromosome];
    parents = reshape(parents, 5,9,[]);
%     parents = parents(:,:,shuffle_ind);
        
    % Crossover
    for j = 1:2:size(parent_bots,2)
        if j== size(parent_bots,2)
            [children(:,:,j), ~] = crossoverMat(parents(:,:,j), parents(:,:,1));
        else
            [children(:,:,j), children(:,:,j+1)] = crossoverMat(parents(:,:,j), parents(:,:,j+1));
        end
    end
    
    %Mutation
    for j = 1:ceil(m*p)
        rand_idx = randi(size(children,3));
        children(rand_idx) = mutate1(genes(rand_idx));
    end
    children_bots = MorphCube(children, [parent_bots.age] + 1);
    
    %Add new random individuals
    rand_bots = MorphCube(rand(5,9,ceil(r*p)));

    %Evaluate
    [child_fits, n_eval_child] = evaluate(sim, children_bots);
    [rand_fits, n_eval_rand] = evaluate(sim, rand_bots);
    n_eval = [n_eval (n_eval(end) + n_eval_child + n_eval_rand)]; %#ok<AGROW>
    
    %Record
    par_layers(:,:,i) = front;
    bots = [parent_bots, children_bots, rand_bots];
    fits = [par_fits, child_fits, rand_fits];
    fit_hist(:,i+1) = fits;
    genes = cat(3, parents, children, reshape([rand_bots.chromosome],5,9,[]));
    divMat(:,:,i+1) = std(genes, 0, 3);
    
    disp(['Completed Gen ' num2str(i)])
end
toc
%%
disp('Done!!');
save('test_run_5');

%%
[M,I] = max(fits);
bot_no = I;

figure;
scatter(par_layers(:,2,end), par_layers(:,1,end), '+')
title('Pareto layer (Last Gen)')

% diversity plot
div = sum(sum(divMat, 2),1);
figure; plot(reshape(div,1,[]), 'm');
title('Diversity Plot');
xlim([1 size(divMat, 3)]);
xlabel('No. of Generations')
ylabel('Sum of Standard Deviation')

% dot plot
figure;
var = reshape(repmat(1:(g + 1), p, 1), [], 1);
scat = scatter(var, reshape(fit_hist(:,1:(g + 1)), [], 1), '.');
scat.MarkerEdgeColor = 'b';
title('Dot Plot')
xlim([1 max(var)])
xlabel('No. of Generations')
ylabel('Fitness: COM displacement (m)')

% learning curve
figure;
plot(n_eval, [0 max(fit_hist)])
title('Learning Curve')

% show the best bot
[M,I] = max(fits);
bot_no = I;

bots(bot_no).plotPDF();
bots(bot_no).plotMaterial();

sim = Simulator(MorphCube(bots(bot_no).chromosome));
figure;
sim.drawRobots;

sim = Simulator();
[frames, K, V, COM, fitness] = sim.simulate_and_plot(MorphCube(bots(bot_no).chromosome));
tic
[fits, n_eval_gen] = sim.evaluate(MorphCube(bots(bot_no).chromosome));
toc

% export to video
myVideo = VideoWriter('BestBot_MorphCube.avi');
myVideo.FrameRate = 25;  % Default 30
myVideo.Quality = 100;    % Default 75
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);

%% multiple robots in one video
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
% 
% %% robot zoo plotting
% bot_zoo_chrom = reshape([bots(10:18).chromosome], 5, 9, []);
% bot_zoo = MorphCube(bot_zoo_chrom);
% 
% figure;
% for i = 1:length(bot_zoo)
%     subplot(3,3,i);
%     sim = Simulator(MorphCube(bot_zoo(i).chromosome));
%     sim.drawRobots()
% end
% 
% bot_zoo(2).plotMaterial;
% 
% example_pdf_bot = MorphCube(bot_zoo(2).chromosome);
% example_pdf_bot.plotPDF()

