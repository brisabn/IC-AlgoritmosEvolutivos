% Configurações iniciais
max_evaluations = 5000;
history = zeros(1, max_evaluations);
learning_rate = 0.01; % Taxa de aprendizado
tolerance = 1e-6; % Tolerância de convergência
x = rand(1, 2); % Ponto inicial aleatório (2D neste exemplo)

% Função objetivo
obj_funct = @(x) funct(x);

for evaluations = 1:max_evaluations
    % Calcule o gradiente numericamente
    gradient = calculate_numeric_gradient(x, obj_funct);
    
    % Atualize a posição usando o Gradiente Descendente
    x = x - learning_rate * gradient;
    
    % Avalie a função no novo ponto
    current_value = obj_funct(x);
    
    % Atualize o histórico
    history(evaluations) = current_value;
    
    % Verifique a convergência
    if evaluations > 1 && abs(history(evaluations) - history(evaluations - 1)) < tolerance
        break; % Convergiu
    end
end

best_fit = history(evaluations);

fprintf('Best value found: %f\n', best_fit);

figure;
plot(history(1:evaluations));
xlabel('Iteration');
ylabel('Best Value');
title('Convergence Plot');

% Função para calcular o gradiente numericamente
function gradient = calculate_numeric_gradient(x, funct)
    epsilon = 1e-6;
    n = length(x);
    gradient = zeros(1, n);
    
    for i = 1:n
        x_plus = x;
        x_minus = x;
        x_plus(i) = x_plus(i) + epsilon;
        x_minus(i) = x_minus(i) - epsilon;
        
        gradient(i) = (funct(x_plus) - funct(x_minus)) / (2 * epsilon);
    end
end
