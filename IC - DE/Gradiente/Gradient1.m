% Função objetivo
obj_funct = @(x) funct(x);

% Pré-alocar history
max_evaluations = 5000;
history = zeros(1, max_evaluations);

% Taxa de aprendizado
learning_rate = 0.01;

% Ponto inicial
x = rand(1, 2); % Inicialize com um ponto aleatório

for evaluations = 1:max_evaluations
    % Avalie a função objetivo no ponto atual
    f = obj_funct(x);
    
    % Calcule o gradiente da função objetivo em relação a x manualmente
    grad = calculate_gradient(obj_funct, x);
    
    % Atualize x usando o Gradiente Descendente
    x = x - learning_rate * grad;
    
    % Registre o valor da função objetivo na iteração atual
    history(evaluations) = f;
    
    % Verifique a convergência
    if evaluations > 1 && abs(history(evaluations) - history(evaluations-1)) < 1e-6
        break; % Convergência detectada
    end
end

best_fit = history(evaluations);

fprintf('Best value found: %f\n', best_fit);

figure;
plot(history(1:evaluations));
xlabel('Iteration');
ylabel('Best Value');
title('Convergence Plot');

function [grad] = calculate_gradient(f, x)
    epsilon = 1e-6;
    grad = zeros(size(x));
    for i = 1:length(x)
        x_plus_epsilon = x;
        x_plus_epsilon(i) = x_plus_epsilon(i) + epsilon;
        grad(i) = (f(x_plus_epsilon) - f(x)) / epsilon;
    end
end
