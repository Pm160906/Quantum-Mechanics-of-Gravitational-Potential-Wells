clc; close all; clear;

%% Step 1: Constants
hbar = 1.054571817e-34;     % Reduced Planck constant [J·s]
G = 6.67430e-11;            % Gravitational constant [m^3/kg·s^2]
m1 = 1.898e27;              % Mass of Jupiter [kg]
m2 = 1.4819e23;             % Mass of Ganymede [kg]

mu = (m1 * m2) / (m1 + m2);         % Reduced mass [kg]
alphas = G * m1 * m2;                % |-alpha| = G * m1 * m2
a_val = hbar^2 / (2 * mu * alphas);  % Bohr-like gravitational radius [m]

fprintf('Gravitational Bohr radius a = %.3e meters\n\n', a_val);

%% Step 2: Symbolic definitions
syms x a v real
assume([a > 0, v > 0])  % Prevent negative radii

psi = sym('psi', [3,1]);       % Unnormalized ψ_n(x)
I_expr = sym('I', [3,1]);      % I(n, L = a, a)

for n = 1:3
    Ln = laguerreL(n-1, 1, v);  % Associated Laguerre polynomial
    R = int(v^2 * exp(-v) * Ln^2, v, 1/n, Inf);  % L = a ⇒ L/(na) = 1/n
    I_expr(n) = simplify(2*n*(n+1) - R);
    
    vsub = x / (n*a);  % Substitution for normalized variable
    prefactor = a^(-3/2) * (n^3 * I_expr(n))^(-1/2);
    psi(n) = simplify(prefactor * (x/a) * exp(-x/(2*n*a)) * laguerreL(n-1, 1, vsub));
end

T = diag(subs(psi, a, a_val));  % Diagonal matrix with unnormalized ψₙ(x)

%% Step 3: Substitute a = a_val and normalize ψₙ(x) over [0, a]
C = sym('C', [3,1]);          % Normalization constants
psi_norm = sym('psi_norm', [3,1]);

for n = 1:3
    psi_sub = simplify(subs(psi(n), a, a_val));
    C(n) = vpa(1 / sqrt(int(psi_sub^2, x, 0, a_val)));  % Compute and store C_n
    psi_norm(n) = simplify(C(n) * psi_sub);             % Normalize the wavefunction
end

%% Step 4: Display results
fprintf('\nDiagonal matrix T (unnormalized ψₙ(x)):\n');
pretty(vpa(T, 4))

fprintf('\nNormalization constants C_n (with L = a):\n');
for n = 1:3
    fprintf('C_%d = %.4e\n', n, C(n));
end

fprintf('\nNormalized wavefunctions ψₙ(x):\n');
digits(4);
for n = 1:3
    fprintf('ψ_%d(x) =\n', n);
    pretty(vpa(psi_norm(n)));
end

%% Step 5: Plot
% Wave function vs x
figure;
fplot(matlabFunction(psi_norm(1)), [0, 20*a_val], 'LineWidth', 1.5)
hold on
fplot(matlabFunction(psi_norm(2)), [0, 20*a_val], 'LineWidth', 1.5)
fplot(matlabFunction(psi_norm(3)), [0, 20*a_val], 'LineWidth', 1.5)

legend('ψ_1(x)', 'ψ_2(x)', 'ψ_3(x)')
xlabel('x (m)')
ylabel('ψ_n(x)')
title('Normalized Gravitational Wavefunctions')
grid on;

% Probability density vs x
% Discretize x
x_vals = linspace(0, 20*a_val, 100);      % 100 spatial points

figure;
for n = 1:3
    % Convert symbolic ψₙ(x) to numeric function
    f_psi = matlabFunction(psi_norm(n), 'Vars', x);
    
    % Evaluate on x_vals
    psi_n_vals = f_psi(x_vals);
    
    % Compute probability density
    prob_density = abs(psi_n_vals).^2;
    
    % Plot
    subplot(3,1,n);
    plot(x_vals, prob_density, 'LineWidth', 1.5);
    xlabel('x (m)');
    ylabel(['|\psi_', num2str(n), '(x)|^2']);
    title(['Probability Density of \psi_', num2str(n), '(x)']);
    grid on;
end

%% Step 6: Data Reduction using SVD
num_x = length(x_vals);
psi_vals = zeros(3, num_x);               % Matrix to hold ψₙ(x) samples

% Evaluate each normalized wavefunction numerically
for n = 1:3
    f = matlabFunction(psi_norm(n), 'Vars', x);  % Convert symbolic to numeric
    psi_vals(n, :) = f(x_vals);                 % Store sampled values
end

% Apply SVD
[U, S, V] = svd(psi_vals, 'econ');   % 'econ' is efficient for non-square

% Extract dominant mode (rank-1 approximation)
rank1 = S(1, 1) * U(:, 1) * V(:, 1)';   % 3×1 * 1×100 ⇒ 3×100 approx matrix

% Visualize weights (U(:,1)) — tells how much ψ₁, ψ₂, ψ₃ contribute
disp('Weights (α_n) for dominant mode (from U(:, 1)):');
disp(U(:, 1));

% Compare original vs approximated ψₙ(x)
figure;
for n = 1:3
    subplot(3,1,n)
    plot(x_vals, psi_vals(n,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(x_vals, rank1(n,:), 'r--', 'LineWidth', 1.5);
    legend(sprintf('Original ψ_%d(x)', n), sprintf('Approx. ψ_%d(x)', n));
    xlabel('x'); ylabel(sprintf('ψ_%d(x)', n));
    title(sprintf('Comparison for ψ_%d(x)', n));
    grid on;
end

%% Step 7: Visualize dominant mode φ(x) and weight αₙ for each ψₙ

% Extract dominant spatial pattern
phi = V(:,1);              % 100×1 ⇒ shape of dominant shared pattern φ(x)
phi = phi / max(abs(phi)); % Normalize for easier viewing

% Extract weights αₙ = s₁ × u₁(n)
alphas = S(1,1) * U(:,1);   % 3×1 vector of contributions to φ(x)

% Display αₙ values
fprintf('\nWeights αₙ (contribution of each ψₙ to φ(x)):\n');
for n = 1:3
    fprintf('α_%d = %.4f\n', n, alphas(n));
end

% Superpose ψₙ(x) using the extracted αₙ values
psi_superposed = alphas(1)*psi_vals(1,:) + alphas(2)*psi_vals(2,:) + alphas(3)*psi_vals(3,:);

% Normalize the superposed wavefunction
norm_super = trapz(x_vals, abs(psi_superposed).^2);
psi_superposed = psi_superposed / sqrt(norm_super);

% Plot it
figure;

plot(x_vals, psi_superposed, 'k-', 'LineWidth', 2);
xlabel('x'); ylabel('\psi_{superposed}(x))');
title('Superposed Wavefunction');
grid on;

%% Step 8: Optimization – Find αₙ minimizing variance of superposed ψ(x)

% Initial guess
alpha0 = ones(3,1)/sqrt(3); % normalized initial guess

% Objective function handle
obj_fun = @(alpha) calcVariance(alpha, psi_vals, x_vals);

% Nonlinear equality constraint (unit norm)
nonlcon = @(alpha) deal([], sum(alpha.^2) - 1);

% Optimization options
opts = optimoptions('fmincon','Display','iter','Algorithm','sqp');

% Run optimization
[alpha_opt, var_min] = fmincon(obj_fun, alpha0, [], [], [], [], [], [], nonlcon, opts);

fprintf('\nOptimized weights αₙ (minimizing spatial variance):\n');
disp(alpha_opt);

% Plot optimized superposed wavefunction
psi_opt = alpha_opt' * psi_vals;

figure;
plot(x_vals, psi_opt, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 2);
xlabel('x'); ylabel('\psi_{opt}(x)');
title('Optimized Superposed Wavefunction (Min Variance)');
grid on;

% Plot probability density before and after optimiztion
figure;
psi_opt_density = abs(psi_opt).^2;
plot(x_vals, abs(psi_superposed).^2, 'b--', 'LineWidth', 2); hold on;
plot(x_vals, psi_opt_density, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('Probability Density');
legend('SVD Superposed |\psi|^2', 'Optimized Superposed |\psi|^2');
title('Comparison of Superposed Wavefunction Densities');
grid on;