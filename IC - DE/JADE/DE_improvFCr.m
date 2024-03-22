% Parâmetros do DE
N = 50;
dim = 2;
minF = 0.2;
maxF = 0.9;

minCr = 0.80;
maxCr = 0.95;

% alcança exatamente 0 as vezes, e parece que muito mais difícil de
% estagnar

% Função objetivo
obj_funct = @(x) funct(x);

% Pré-alocar history
max_evaluations = 400;
history = zeros(1, max_evaluations);

historyF = zeros(1, max_evaluations);
historyCr = zeros(1, max_evaluations);

% Best inicial
best_fit = inf;
best_sol = zeros(dim, 1);

% Build Population
lower_bound = [-20, -20]';
upper_bound = [20, 20]';
population = lower_bound + (upper_bound - lower_bound) .* rand(dim, N);

% Evaluate all fitness functions of the initial population
population_fitness = zeros(N, 1);
for i = 1:N
    xi = population(:, i);
    population_fitness(i) = obj_funct(xi);
end

% DE loop
for evaluations = 1:max_evaluations
    population2 = zeros(dim, N);
    fitness2 = zeros(N, 1);
    % archives
    archiveF = zeros(1, N);    
    archiveImprov = zeros(1, N);
    archiveCr = zeros(1, N);

    for i = 1:N % Para cada indivíduo da população

        F = rand * (maxF - minF) + minF;
        Cr = rand * (maxCr - minCr) + minCr;

        idx = randperm(N - 1, 3);

        candidates = 1:N;
        candidates(i) = [];

        i2 = candidates(idx(1));
        i3 = candidates(idx(2));
        i4 = candidates(idx(3));

        x1 = population(:, i);
        x2 = population(:, i2);
        x3 = population(:, i3);
        x4 = population(:, i4);

        A = crossover(Cr, dim);
        nA = 1 - A;

        trial = nA.*x1 + A.*(x2 + F*(x3 - x4));
        trial_fitness = obj_funct(trial);

        fitness_old = population_fitness(i); 
        fitness_new = trial_fitness;

        if fitness_new < fitness_old % melhorou
            population2(:, i) = trial;
            fitness2(i) = fitness_new;
            % arquivar melhoria
            improv_step = population_fitness(i) - trial_fitness; % norm(x1 - trial);
            archiveImprov(i) = improv_step;
        else 
            population2(:, i) = x1;
            fitness2(i) = fitness_old;
        end

        archiveF(i) = F; % Salvar o F do indivíduo
        archiveCr(i) = Cr; % Salvar o Cr do indivíduo

        if fitness_new < best_fit
            best_sol = trial;
            best_fit = fitness_new;
        end

    end % Fim da população e encontrou best fit

    history(evaluations) = best_fit;

    % média dos 50% melhores F que deram melhor passo
    new = find(archiveImprov > 0);
    if(~isempty(new))
        [M, I] = sort(archiveImprov(new));

        bestF = archiveF(I(ceil(length(new)/2):length(new)));
        mediaF = mean(bestF);

        bestCr = archiveCr(I(ceil(length(new)/2):length(new)));
        mediaCr = mean(bestCr);
    end 

    % Armazenar o valor de F no historyF
    historyF(evaluations) = mediaF; % pega a última média
    historyCr(evaluations) = mediaCr; % pega a última média

    % novos limites de F
    minF_new = max(0.0001, mediaF - 0.01);
    maxF_new = mediaF + 0.01;

    minF = minF_new;
    maxF = maxF_new;

    % novos limites de Cr
    minCr_new = max(0.0001, mediaCr - 0.01);
    maxCr_new = mediaCr + 0.01;

    minCr = minCr_new;
    maxCr = maxCr_new;

    % Update population and fitness for the next iteration
    population = population2;
    population_fitness = fitness2;

end

fprintf('Best value found: %f at position: (%.20f, %.20f)\n', best_fit, best_sol(1), best_sol(2));
% Verificar se o melhor valor encontrado é exatamente 0
if best_fit == 0
    fprintf('Exatamente 0\n');
end

% Plotar histórico de convergência
% figure;
semilogy(history(1:evaluations));
xlabel('Iteration');
ylabel('Best Value');
title('Convergence Plot');


% Plotar histórico de F
% figure;
% plot(historyF(1:evaluations));
% xlabel('Iteration');
% ylabel('Mean F Value');
% title('Mean F Value Evolution Plot');

% % Filtrar os valores diferentes de 0 para minF e maxF
% nonZero_minF_history = minF_history(minF_history ~= 0);
% nonZero_maxF_history = maxF_history(maxF_history ~= 0);
% 
% figure;
% plot(1:length(nonZero_minF_history), nonZero_minF_history, 'rx', 'MarkerSize', 4); % Pontos vermelhos para minF
% hold on;
% plot(1:length(nonZero_maxF_history), nonZero_maxF_history, 'bx', 'MarkerSize', 4); % Pontos azuis para maxF
% hold off;
% xlabel('Iteration');
% ylabel('Best Value');
% title('Convergence Plot');
% legend('minF', 'maxF');

function A = crossover(Cr, dim)
    A = zeros(dim, 1);
    for d = 1:dim
        if rand() < Cr
            A(d) = 1;  % Trocar
        else
            A(d) = 0;  % Não trocar
        end
    end
end