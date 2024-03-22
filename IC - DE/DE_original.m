% aquele fake que dava 0.2 num e 0.9 noutro sempre
N = 50;
dim = 2;
F = 0.5;  % Valor inicial de F
Cr = 0.5;  % Valor inicial de Cr

obj_funct = @(x) funct(x);

max_evaluations = 5000;
history = zeros(1, max_evaluations);

best_fit = inf;
best_sol = zeros(dim, 1);

lower_bound = [-20, -20]';
upper_bound = [20, 20]';

population = lower_bound + (upper_bound - lower_bound) .* rand(dim, N);

population_fitness = zeros(N, 1);
for i = 1:N
    xi = population(:, i);
    population_fitness(i) = obj_funct(xi);
end

evaluations = 0;

while evaluations < max_evaluations

    population2 = zeros(dim, N);
    fitness2 = zeros(N, 1);

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

        trial = nA .* x1 + A .* (x2 + F * (x3 - x4));
        trial_fitness = obj_funct(trial);
        
        evaluations = evaluations + 1;

        fitness_old = population_fitness(i);
        fitness_new = trial_fitness;

        if fitness_new < fitness_old
            population2(:, i) = trial;
            fitness2(i) = fitness_new;
        else
            population2(:, i) = x1;
            fitness2(i) = fitness_old;
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

    % % Atualize F e Cr com base na convergência
    % if evaluations > 1
    %     % Verifique se a convergência melhorou
    %     if history(evaluations) < history(evaluations - 1)
    %         F = min(0.9, F + 0.1);  % Aumenta F
    %         Cr = max(0.2, Cr - 0.1);  % Diminui Cr
    %     else
    %         F = max(0.2, F - 0.1);  % Diminui F
    %         Cr = min(0.9, Cr + 0.1);  % Aumenta Cr
    %     end
    % end
end

fprintf('Best value found: %f\n', best_fit);
figure;
plot(history(1:evaluations));
xlabel('Iteration');
ylabel('Best Value');
title('Convergence Plot');

function A = crossover(Cr, dim)
    A = zeros(dim, 1);
    for d = 1:dim
        if rand() < Cr
            A(d) = 1;
        else
            A(d) = 0;
        end
    end
end
