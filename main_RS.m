% Random Search
clc
clear
close all
%%
p = 5; % Population size
g = 2; % number of generations

run_time_per_robot = 7; % s
disp(['Estimated run time: ' num2str(floor(run_time_per_robot*p*(g + 1)/3600/2)) ' hr. ' ...
    num2str(round(mod(run_time_per_robot*p*(g + 1)/60/2, 60))) ' mins'])

%%
n_eval = [];
genes= rand(5,9,p);
bots = MorphCube(genes);
sim = Simulator();

divMat = zeros(5,9,g+1);
divMat(:,:,1) = std(genes, 0, 3);

tic
[fits, n_eval_gen] = sim.evaluate(bots);
toc

n_eval = [n_eval n_eval_gen];

fit_hist = [];
fit_hist = [fit_hist, fits'];
[M,I] = max(fits);
best_fit = [M; n_eval];
best_bot = bots(I);

%% start the EA
tic

cont_inc = 0; % specify how many more gens to run 

for i = (1 : g) + cont_inc
    
    %Add new random individuals
    genes = rand(5,9,p);
    bots = MorphCube(genes);

    %Evaluate
    [fits, n_eval_gen] = evaluate(sim, bots);
    n_eval = [n_eval (n_eval(end) + n_eval_gen)]; %#ok<AGROW>
    
      
    % grow the fit_hist 
    fit_hist = [fit_hist, fits']; %#ok<AGROW>
        
    if i <= g
       divMat(:,:,i+1) = std(genes, 0, 3);
    else
       divMat = cat(3, divMat, std(genes, 0, 3)); 
    end
    
    % Save best bot
    [M,I] = max(fits);
    if M > best_fit(1,end)
        best_fit = [best_fit [M; n_eval(end)]];
        best_bot = bots(I);
    end
    
    disp(['Completed Gen ' num2str(i)])
end
toc
disp('Done!!');
%%
save('test_run_9');
%%

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
plot(n_eval, max(fit_hist));
title('Best per gen')

figure;
plot(best_fit(2,:), best_fit(1,:));
title('Learning Curve')

% show the best bot

div = sum(sum(divMat, 2),1);
figure; plot(reshape(div,1,[]));

best_bot.plotPDF();
best_bot.plotMaterial();

sim = Simulator(MorphCube(best_bot.chromosome));
figure;
sim.drawRobots;

sim = Simulator();
[frames, K, V, COM, fitness] = sim.simulate_and_plot(MorphCube(best_bot.chromosome));
tic
[fit, n_eval_gen] = sim.evaluate(MorphCube(best_bot.chromosome));
toc

% export to video
myVideo = VideoWriter('BestBot_RandSearch.avi');
myVideo.FrameRate = 25;  % Default 30
myVideo.Quality = 100;    % Default 75
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);

