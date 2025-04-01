
clear; clc;

%% Parameters
m = 2; % Mass of the hexadrone (kg)
g = 9.81; % Gravity (m/s^2)
c_f = 9.9 * 10^-4; % Lift coefficient
c_d = [-1, +1, -1, +1, -1, +1] * 1.9 * 10^-5; % Drag coefficient

% Initial guesses for parameters
l_init = [0.35, 0.35, 0.35, 0.35, 0.35, 0.35];

gamma_init = deg2rad([0, 60, 120, 180, 240, 300]);
alpha_init = deg2rad([-45, -30, -15, -45, -30, -15]);
beta_init = deg2rad([-45, -30, -15, -45, -30, -15]);
x0 = [l_init, gamma_init, alpha_init, beta_init];

% Lower and Upper Bounds
lb = [0.2 * ones(1,6), deg2rad(0) * ones(1,6), deg2rad(-90) * ones(1,12)];
ub = [1 * ones(1,6), deg2rad(360) * ones(1,6), deg2rad(90) * ones(1,12)];

% Ensure x0 is within bounds
x0 = max(min(x0, ub), lb);

%% Define Optimization Problem
objective = @(x) power_consumption(x, m, g, c_f, c_d);
constraints = @(x) hover_constraints(x, m, g, c_f, c_d);

% Optimization Options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxIterations', 50000, 'MaxFunctionEvaluations', 1000000);

% Run Optimization
[x_opt, fval, exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, constraints, options);

% Display Results
if exitflag > 0
    disp('✅ Optimization Converged Successfully');
    
    % Extract parameters
    l_opt = x_opt(1:6);
    gamma_opt = x_opt(7:12);
    alpha_opt = x_opt(13:18);
    beta_opt = x_opt(19:24);
    
    % Display results in structured way
    fprintf('\n=== OPTIMIZED PARAMETERS ===\n');
    fprintf('Optimal Power Consumption: %.4f\n', fval);
    
    % Display arm lengths
    fprintf('\n--- Arm Lengths (m) ---\n');
    for i = 1:6
        fprintf('Arm %d: %.4f\n', i, l_opt(i));
    end
    
    % Display angles in degrees (more intuitive than radians)
    fprintf('\n--- Gamma Angles (deg) ---\n');
    for i = 1:6
        fprintf('Gamma %d: %.2f\n', i, rad2deg(gamma_opt(i)));
    end
    
    fprintf('\n--- Alpha Angles (deg) ---\n');
    for i = 1:6
        fprintf('Alpha %d: %.2f\n', i, rad2deg(alpha_opt(i)));
    end
    
    fprintf('\n--- Beta Angles (deg) ---\n');
    for i = 1:6
        fprintf('Beta %d: %.2f\n', i, rad2deg(beta_opt(i)));
    end
    
    % Calculate rotor speeds at optimal point
    [F, efficiency] = getRotorSpeed(alpha_opt, beta_opt, gamma_opt, l_opt);
    w_opt = sqrt(max(abs(pinv(F) * [0; 0; m * g; 0; 0; 0]), 0));
    
    fprintf('\n--- Rotor Speeds ---\n');
    for i = 1:6
        fprintf('Rotor %d: %.2f\n', i, w_opt(i));
    end
    
    fprintf('\nEfficiency: %.4f\n', efficiency);
    
    % Calculate constraint values
    [c, ~] = hover_constraints(x_opt, m, g, c_f, c_d);
    c1 = c(1:6);  % Lower bound constraints
    c2 = c(7:12); % Upper bound constraints
    
    fprintf('\n--- Constraint Violations ---\n');
    fprintf('Min speed violation (should be <= 0): %.4f\n', max(c1));
    fprintf('Max speed violation (should be <= 0): %.4f\n', max(c2));
    
    % Optional: Visualize the hexadrone configuration
    try
        visualize_hexadrone(l_opt, gamma_opt, alpha_opt, beta_opt);
    catch
        fprintf('\nVisualization requires additional function.\n');
    end
else
    disp('❌ Optimization Failed');
end


%% Objective Function: optimitzig for the power consumption
function P = power_consumption(x, m, g, c_f, c_d)
    l = x(1:6);
    gamma = x(7:12);
    alpha = x(13:18);
    beta = x(19:24);

    % Compute Rotor Speeds
    [F, ~] = getRotorSpeed(alpha, beta, gamma, l);
    % Dealung with Singular Matrix Issue
    if rank(F) < 6 || cond(F) > 1e10 % Avoid near-singular matrices
        P = inf; % Penalize invalid solutions
        return;
    end
    % Use Pseudo-Inverse Instead of inv(F) in case of siongularity
    w = sqrt(max(abs(pinv(F) * [0; 0; m * g; 0; 0; 0]), 0));
    % Check for NaN values
    if any(isnan(w)) || any(isinf(w))
        P = inf; % Penalize bad solutions
        return;
    end
    P = sum(w.^2);
end
%% Constraint Function with Minimum Gamma Spacing;
%   in this version of the fucntion, a minimum spacing for the angle gamma
%   is introduced, as without the spacing the optimizer converges to
%   unfeasable solutions in which most of the rotors are aligned on a
%   single axis.
%   the without spacing function was kept ( scroll down) for comparision
%   reasons.

function [c, ceq] = hover_constraints(x, m, g, c_f, c_d)
    l = x(1:6);
    gamma = x(7:12);
    alpha = x(13:18);
    beta = x(19:24);
    
    % Initialize constraints
    ceq = []; % No equality constraints
    [F, ~] = getRotorSpeed(alpha, beta, gamma, l);
    
    % Check for rank/condition issues
    if rank(F) < 6 || cond(F) > 1e10
        % Need to return correct number of constraints even in this case
        c = ones(27, 1) * 1000; % 12 original + 15 new gamma spacing constraints
        return;
    end
    
    % Compute rotor speeds
    w = sqrt(max(abs(pinv(F) * [0; 0; m * g; 0; 0; 0]), 0));
    
    % Ensure rotor speeds are within bounds (25 ≤ w ≤ 100)
    c1 = 25 - w; % 6 constraints
    c2 = w - 100; % 6 constraints
    
    % Add minimum spacing constraints between all pairs of gamma values
    % 20 degrees in radians
    min_spacing = deg2rad(23);
    gamma_constraints = [];
    
    % Check all pairs of gamma values
    for i = 1:5
        for j = (i+1):6
            % Calculate absolute angular difference, considering periodicity
            diff = abs(angdiff(gamma(i), gamma(j)));
            % The constraint is violated if diff < min_spacing
            gamma_constraints = [gamma_constraints; min_spacing - diff];
        end
    end
    
    % Combine all constraints
    c = [c1; c2; gamma_constraints];
end


% Comstraint funciton without the spacing between the gamma angles; unused 

function [c, ceq] = hover_constraints_without_spacing(x, m, g, c_f, c_d)
    l = x(1:6);
    gamma = x(7:12);
    alpha = x(13:18);
    beta = x(19:24);
    % Initialize constraints
    ceq = []; % No equality constraints
    [F, ~] = getRotorSpeed(alpha, beta, gamma, l);
    % Check for rank/condition issues
    if rank(F) < 6 || cond(F) > 1e10
        % Need to return 12 elements even in this case
        c = zeros(12, 1);
        c(:) = 1000; % Large penalty value
        return;
    end
    % Compute rotor speeds
    w = sqrt(max(abs(pinv(F) * [0; 0; m * g; 0; 0; 0]), 0));
    % Ensure rotor speeds are within bounds (25 ≤ w ≤ 100)
    c1 = 25 - w; % 6 constraints
    c2 = w - 100; % 6 constraints
    c = [c1; c2]; % Total of 12 constraints
end

%% helper function if angdiff is not available
function d = angdiff(a1, a2)
    % Returns the smallest angle difference between two angles
    % considering the periodicity of angles
    d = mod(a1 - a2 + pi, 2*pi) - pi;
    d = abs(d);
end
%% Rotor Speed Calculation Function
function [F, efficiency] = getRotorSpeed(alpha, beta, gamma, l)
    c_f = 9.9 * 10^-4;
    c_d = [-1, +1, -1, +1, -1, +1] * 1.9 * 10^-5;
    F1 = zeros(3, 6);
    F2 = zeros(3, 6);
    for i = 1:6
        Rb_S = Rz(gamma(i)) * Ry(beta(i)) * Rx(alpha(i));
        F1(:, i) = Rb_S * [0; 0; c_f];
        F2(:, i) = cross(Rz(gamma(i)) * [l(i); 0; 0], Rb_S * [0; 0; c_f]) + Rb_S * [0; 0; c_d(i)];
    end
    efficiency = sum(F1(3, :)) / (6 * c_f);
    F = [F1; F2];
end
%% Rotation Matrices
function Rx = Rx(alpha)
    Rx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
end
function Ry = Ry(beta)
    Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
end
function Rz = Rz(gamma)
    Rz = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
end
%% Optional Visualization Function
function visualize_hexadrone(l, gamma, alpha, beta)
    figure('Name', 'Hexadrone Configuration', 'Position', [100, 100, 800, 600]);
    
    % Origin (drone center)
    plot3(0, 0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    hold on;
    
    % Color map for rotors
    colors = {'r', 'g', 'b', 'c', 'm', 'y'};
    
    % Plot each arm and rotor
    for i = 1:6
        % Arm endpoint in global frame
        arm_end = Rz(gamma(i)) * [l(i); 0; 0];
        
        % Draw arm
        plot3([0, arm_end(1)], [0, arm_end(2)], [0, arm_end(3)], '-', 'Color', colors{i}, 'LineWidth', 2);
        
        % Rotor orientation
        rotor_normal = Rz(gamma(i)) * Ry(beta(i)) * Rx(alpha(i)) * [0; 0; 1];
        rotor_normal = rotor_normal / norm(rotor_normal) * 0.1; % Scale for visualization
        
        % Draw rotor normal vector
        quiver3(arm_end(1), arm_end(2), arm_end(3), ...
                rotor_normal(1), rotor_normal(2), rotor_normal(3), ...
                'Color', colors{i}, 'LineWidth', 2);
        
        % Plot rotor position
        plot3(arm_end(1), arm_end(2), arm_end(3), 'o', 'Color', colors{i}, ...
              'MarkerSize', 8, 'MarkerFaceColor', colors{i});
        
        % Add text label
        text(arm_end(1), arm_end(2), arm_end(3), ['  Rotor ', num2str(i)], 'Color', colors{i});
    end
    
    % Add coordinate axes
    axis_length = 1.5 * max(l);
    quiver3(0, 0, 0, axis_length, 0, 0, 'k', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, axis_length, 0, 'k', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 0, axis_length, 'k', 'LineWidth', 2);
    text(axis_length, 0, 0, 'X', 'FontSize', 12);
    text(0, axis_length, 0, 'Y', 'FontSize', 12);
    text(0, 0, axis_length, 'Z', 'FontSize', 12);
    
    grid on;
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Optimized Hexadrone Configuration');
    view(30, 30);
end

