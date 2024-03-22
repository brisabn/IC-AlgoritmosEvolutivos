% Parâmetros do DE
N = 50;
dim = 2;
F = 0.9;
Cr = 0.8;

% Função objetivo
obj_funct = @(x) funct(x);

% Pré-alocar history
max_evaluations = 200;
history = zeros(1, max_evaluations);

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

    % Encontre o índice do melhor indivíduo
    [~, best_idx] = min(population_fitness);

    for i = 1:N
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

        if fitness_new < fitness_old
            population2(:, i) = trial;
            fitness2(i) = fitness_new;
        else
            population2(:, i) = x1;
            fitness2(i) = fitness_old;
        end

        % Inclua o melhor indivíduo da geração anterior na nova população
        if fitness_new < population_fitness(best_idx)
            population2(:, best_idx) = population(:, best_idx);  % Manter o melhor indivíduo anterior
            fitness2(best_idx) = population_fitness(best_idx);
        else
            population2(:, best_idx) = best_sol;  % Inserir o melhor indivíduo anterior
            fitness2(best_idx) = best_fit;
        end

        if fitness_new < best_fit
            best_sol = trial;
            best_fit = fitness_new;
        end

        history(evaluations) = best_fit;
    end

    % Update population and fitness for the next iteration
    population = population2;
    population_fitness = fitness2;
end

fprintf('Best value found: %f at position: (%.20f, %.20f)\n', best_fit, best_sol(1), best_sol(2));
figure;
semilogy(history(1:evaluations));
xlabel('Iteration');
ylabel('Best Value');
title('Convergence Plot');

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
