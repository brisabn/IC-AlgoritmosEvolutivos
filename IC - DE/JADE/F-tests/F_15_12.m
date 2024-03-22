% Parâmetros do DE
N = 50;
dim = 2;
minF = 0.2;
maxF = 0.9;
Cr = 0.8;

% Função objetivo
obj_funct = @(x) funct(x);

% Pré-alocar history
max_evaluations = 200;
history = zeros(1, max_evaluations);
minF_history = zeros(1, max_evaluations);
maxF_history = zeros(1, max_evaluations);
historyF = zeros(1, max_evaluations);

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
    
    for i = 1:N % Para cada indivíduo da população

        F = rand * (maxF - minF) + minF;

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
            % arquiva o improvement melhor
            improv_step = population_fitness(i) - trial_fitness; % norm(x1 - trial);
            archiveImprov(i) = improv_step;
        else 
            population2(:, i) = x1;
            fitness2(i) = fitness_old;
        end

        archiveF(i) = F; % Salvar o F do indivíduo

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
        media = mean(bestF);
    end 

    % Armazenar o valor de F no historyF
    historyF(evaluations) = media; % pega a última média

    % novos limites
    minF_new = max(0.0001, media - 0.01);
    maxF_new = media + 0.01;
    
    % interval_size = max(0.2, maxF - minF);
    % minF_new = max(0.2, media - interval_size/2);
    % maxF_new = media + interval_size/2;

    % atualizar os limites de F
    minF = minF_new;
    maxF = maxF_new;
    minF_history(evaluations) = minF;
    maxF_history(evaluations) = maxF;

    % Update population and fitness for the next iteration
    population = population2;
    population_fitness = fitness2;

end

fprintf('Best value found: %f\n', best_fit);

% Plotar histórico de convergência
figure;
semilogy(history(1:evaluations));
xlabel('Iteration');
ylabel('Best Value');
title('Convergence Plot');

% Plotar histórico de F
figure;
plot(historyF(1:evaluations));
xlabel('Iteration');
ylabel('Mean F Value');
title('Mean F Value Evolution Plot');

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