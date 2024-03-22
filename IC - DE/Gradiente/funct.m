function y = funct(x)

n = length(x);
fun = 7;
Q = matpi(n);

if fun==1  % Quadratic function
    phi = 1e-3;
    v1 = ones(n,1)/sqrt(n);
    v2 = null(v1');
    V = [v1 v2];
    d = ones(1,n);
    d(1) = phi;
    D = diag(d);
    Q = V'*D*V/phi;
    y = x'*Q*x;
elseif fun==2  % Sphere function
    y = x'*x;
elseif fun==3 % Rosenbrock function
    y = 0;
    for i=1:n-1
        y = y + 100*(x(i+1) - x(i)^2)^2 + (x(i) - 1)^2;
    end
elseif fun==4  % rotated Rastrigin function
    z = Q*x;
    y = 0;
    for i=1:n
        y = y + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
    end
elseif fun==5  % rotated Ackley function
    z = Q*x;
    cz = cos(2*pi*z);
    y = -20*exp(-0.2*sqrt(z'*z/n)) - exp(sum(cz)/n) + 20 + exp(1);
elseif fun==6  % Schwefel function
    y = 0;
    for i=1:n
        y = y + (sum(x(1:i)))^2;
    end
elseif fun==7  % Your new function
    y = 0;
    for k = -10:10
        y = y + exp(-(x(1)-x(2))^2 - 2*x(1)^2)*cos(x(2))*sin(2*x(2));
    end
end


% ================================================================
function M = matpi(n)
% matriz n x n contendo os digitos das potencias inteiras de pi
% M é uma matriz construída de maneira determinística

npi = ceil(n^2/15);
a = [];
for i=1:npi
    a = [a num2str(pi^i,17)];
end
posp = find(a=='.');
a(posp) = [];
M = zeros(n);
for i=1:n
    for j=1:n
        ind = (i-1)*n + j;
        M(i,j) = str2num(a(ind));
    end
end
M = M'*M;
M = M/norm(M);

