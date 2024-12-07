clear; clc; close all;

% Initial conditions
masses = [1, 1, 1]; % Example masses
y0 = [
    -0.97000436;  0.24308753;  % Initial position of body 1
     0.46620368;  0.43236573;  % Initial velocity of body 1
     0.97000436; -0.24308753;  % Initial position of body 2
     0.46620368;  0.43236573;  % Initial velocity of body 2
     0;          0;           % Initial position of body 3
    -0.93240737; -0.86473146;  % Initial velocity of body 3
];
tspan = [0 100]; % Time range for the simulation
dt = 0.01;
num_steps = 1000; % Number of prediction steps
degree = 3;      % Degree of polynomial basis
figTitle = 'Three-Body Problem Predicted Trajectories Using EDMD';
% Initialize array to store computation times
computation_times = zeros(1, 4); % 4 methods: Euler, ODE45, RKF4(5), Leapfrog

plotInitialConditions(y0);

% Solve the system using built in matlab ODE45 (Runge-Kutta) solver
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
tic;
[timeRK, yRK] = ode45(@(t, y) threeBodyODE(t, y, masses), tspan, y0);
computation_times(2) = toc;
plotPhaseSpace(yRK, "Runge Kutta Method");

% Solve using Euler method
tic; % Start timer
[timeE, yE] = eulerThreeBody(@threeBodyODE, masses, y0, tspan, dt);
computation_times(1) = toc; % Stop timer
plotPhaseSpace(yE, "Euler Method");

% Solve using RKF4(5)
tol = 1e-6;
tic;
[timeRKF, yRKF] = rkf45ThreeBody(@threeBodyODE, masses, y0, tspan, tol);
computation_times(3) = toc;
plotPhaseSpace(yRKF, "RKF4(5) Method");

% Solve using Leapfrog
dtLF = 0.001;
tic;
[timeLF, yLF] = leapfrogThreeBody(@threeBodyODE, masses, y0, tspan, dtLF);
computation_times(4) = toc;
plotPhaseSpace(yLF, "Leapfrog");

% Display results
methods = {'Euler', 'ODE45', 'RKF4(5)', 'Leapfrog'};
disp('Computation Times for Each Method:');
for i = 1:length(methods)
    fprintf('%s: %.6f seconds\n', methods{i}, computation_times(i));
end

% Plot computation times
figure;
bar(categorical(methods), computation_times);
ylabel('Computation Time (seconds)');
title('Computation Time Comparison Across Methods');
grid on;

% Combine into cell arrays
t_values = {timeE, timeRK, timeRKF, timeLF};
y_values = {yE, yRK, yRKF, yLF};

% Visualise the ODE45 solution (e.g., positions of the three bodies)
r1 = yRK(:, 1:2);
r2 = yRK(:, 5:6);
r3 = yRK(:, 9:10);
figure;
plot(r1(:, 1), r1(:, 2), 'r', r2(:, 1), r2(:, 2), 'g', r3(:, 1), r3(:, 2), 'b');
xlabel('x'); ylabel('y'); legend('Body 1', 'Body 2', 'Body 3');
title('Three-Body Problem Trajectories (Runge-Kutta)');

dtRK = diff(timeRK); % Compute the adaptive time steps
figure;
plot(timeRK(1:end-1), dtRK);
xlabel('Time');
ylabel('Time Step (\Delta t)');
title('Adaptive Time Step Sizes (Runge Kutta)');

% Visualise the Euler method solution (e.g., positions of the three bodies)
r1e = yE(:, 1:2);
r2e = yE(:, 5:6);
r3e = yE(:, 9:10);
figure;
plot(r1e(:, 1), r1e(:, 2), 'r', r2e(:, 1), r2e(:, 2), 'g', r3e(:, 1), r3e(:, 2), 'b');
xlabel('x'); ylabel('y'); legend('Body 1', 'Body 2', 'Body 3');
title('Three-Body Problem Trajectories (Euler Method)');

% Visualise the RKF4(5) method solution (e.g., positions of the three bodies)
r1RKF = yRKF(:, 1:2);
r2RKF = yRKF(:, 5:6);
r3RKF = yRKF(:, 9:10);
figure;
plot(r1RKF(:, 1), r1RKF(:, 2), 'r', r2RKF(:, 1), r2RKF(:, 2), 'g', r3RKF(:, 1), r3RKF(:, 2), 'b');
xlabel('x'); ylabel('y'); legend('Body 1', 'Body 2', 'Body 3');
title('Three-Body Problem Trajectories (RKF4(5))');

dtRKF = diff(timeRKF); % Compute the adaptive time steps
figure;
plot(timeRKF(1:end-1), dtRKF);
xlabel('Time');
ylabel('Time Step (\Delta t)');
title('Adaptive Time Step Sizes (RKF4(5))');

% Visualise the RKF4(5) method solution (e.g., positions of the three bodies)
r1LF = yLF(:, 1:2);
r2LF = yLF(:, 5:6);
r3LF = yLF(:, 9:10);
figure;
plot(r1LF(:, 1), r1LF(:, 2), 'r', r2LF(:, 1), r2LF(:, 2), 'g', r3LF(:, 1), r3LF(:, 2), 'b');
xlabel('x'); ylabel('y'); legend('Body 1', 'Body 2', 'Body 3');
title('Three-Body Problem Trajectories (Leapfrog)');

%Truncation error evaluation (Euler)
dt_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005];
evaluateTruncationError(@threeBodyODE, masses, timeLF, yLF, y0, tspan, dt_values);

periodic_differences = zeros(size(dt_values));

for i = 1:length(dt_values)
    dt = dt_values(i);
    [~, y] = eulerThreeBody(@threeBodyODE, masses, y0, tspan, dt);
    % Distance between initial and final states
    periodic_differences(i) = norm(y(end, :) - y(1, :));
end

% Plot periodic differences
figure;
plot(dt_values, periodic_differences, '-o', 'LineWidth', 2);
xlabel('Step size (\Delta t)');
ylabel('Periodic Difference (\|y(T) - y(0)\|)');
title('Periodicity Check for Different \Delta t');
grid on;


% Compute energy for all methods
energyEuler = computeEnergy(yE, masses);
energyRK = computeEnergy(yRK, masses);
energyRKF = computeEnergy(yRKF, masses);
energyLF = computeEnergy(yLF, masses);

% Plot energy conservation for all methods
figure;
plot(timeE, energyEuler, 'k', 'LineWidth', 1.5); hold on;
plot(timeRK, energyRK, 'b', 'LineWidth', 1.5);
plot(timeRKF, energyRKF, 'r--', 'LineWidth', 1.5);
plot(timeLF, energyLF, 'g-.', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Total Energy');
title('Energy Conservation Comparison for All Methods');
legend('Euler', 'ODE45', 'RKF4(5)', 'Leapfrog', 'Location', 'best');
grid on;

% Plot energy conservation for leapfrog and rkf4(5)
figure;
plot(timeLF, energyLF, 'g-.', 'LineWidth', 1.5); hold on;
plot(timeRKF, energyRKF, 'r--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Total Energy');
title('Energy Conservation Comparison for Leapfrog and RKF4(5)');
legend('Leapfrog', 'RKF4(5)');
grid on;

% Plot energy conservation for leapfrog
figure;
plot(timeLF, energyLF, 'g-.', 'LineWidth', 1.5); hold on;
xlabel('Time');
ylabel('Total Energy');
title('Energy Conservation Comparison for Leapfrog');
legend('Leapfrog');
grid on;

methods = {'Euler', 'ODE45', 'RKF4(5)', 'Leapfrog'};
% Call the periodicity comparison function
plotPeriodicityComparison(methods, t_values, y_values, y0);



% Compute the distance to the initial state
distances = vecnorm(yRK - y0', 2, 2); % L2 norm of the difference

% Plot the periodicity
figure;
plot(timeRK, distances, 'Color', "r", 'LineWidth', 1.5, ...
'DisplayName', "Runge Kutte Periodicity");
grid on;

% Step 1: Extract the Envelope
[~, locs] = findpeaks(-distances); % Find local maxima (upper envelope)
x_env = timeRK(locs);          % x-values of the envelope
y_env = abs(distances(locs));     % y-values of the envelope (absolute for symmetry)

threshold = 2; % Adjust as needed
filtered_indices = y_env < threshold; % Logical indices where troughs are below threshold
x_filtered = x_env(filtered_indices); % Filtered x-values
y_filtered = y_env(filtered_indices); % Filtered y-values

% Step 2: Log-transform the Envelope
log_y_env = log(y_filtered); % Apply natural logarithm

% Step 3: Fit a Linear Model to Log-Transformed Data
p = polyfit(x_filtered, log_y_env, 1); % Linear fit: log(y) = log(A) - k*x
k = -p(1); % Decay rate (negative slope)
log_a = p(2); % log(A)
a = exp(log_a); % Initial amplitude

% Step 4: Visualize Results
figure;

% Original data
subplot(2, 1, 1);
plot(timeRK, distances, 'b-', 'DisplayName', 'Original Data'); hold on;
plot(x_filtered, y_filtered, 'ro', 'MarkerSize', 8, 'DisplayName', 'Envelope Points');
plot(timeRK, a * exp(-k * timeRK), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Fitted Exponential Envelope');
xlabel('timeRK');
ylabel('distances');
legend('show');
title('Original Data and Fitted Exponential Envelope');
grid on;

% Log-transformed envelope
subplot(2, 1, 2);
plot(x_filtered, log_y_env, 'ro-', 'DisplayName', 'Log of Envelope'); hold on;
plot(x_filtered, polyval(p, x_filtered), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit');
xlabel('timeRK');
ylabel('log(distances)');
legend('show');
title('Log-Transformed Envelope and Linear Fit');
grid on;

% Step 5: Display Fit Results
fprintf('Exponential Decay Fit:\n');
fprintf('Initial Amplitude (A): %.4f\n', a);
fprintf('Decay Rate (k): %.4f\n', k);


%Animate Trajectories
%animateThreeBody(timeRK, yRK);
%animateThreeBody(timeE, yE);
%animateThreeBody(timeRKF, yRKF);
animateThreeBody(timeLF, yLF);

%Visualise Initial Conditions
function plotInitialConditions(y0)
    % plotInitialConditions: Plots the initial conditions of the 3-body problem
    % with velocity vectors.
    %
    % Inputs:
    %   y0: Initial conditions in the format:
    %       [x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3]

    % Reshape y0 into a 2D array for easier indexing
    y0 = reshape(y0, 4, 3)'; % Each row corresponds to a body [x, y, vx, vy]
    
    % Extract positions and velocities
    positions = y0(:, 1:2);  % Positions of the three bodies
    velocities = y0(:, 3:4); % Velocities of the three bodies
    
    % Create the plot
    figure;
    hold on;
    grid on;
    
    % Plot positions and velocity vectors
    for i = 1:3
        % Plot body positions
        plot(positions(i, 1), positions(i, 2), 'o', 'MarkerSize', 8, ...
             'DisplayName', ['Body ', num2str(i)]);
        
        % Plot velocity vectors
        quiver(positions(i, 1), positions(i, 2), velocities(i, 1), velocities(i, 2), ...
            'AutoScale', 'off', 'MaxHeadSize', 0.5, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', ['Velocity ', num2str(i)]);
    end
    
    % Formatting
    xlabel('x');
    ylabel('y');
    title('Initial Conditions of the 3 Bodies with Velocity Vectors');
    legend;
    axis equal;
    hold off;
end

%Phase Space
function plotPhaseSpace(y, figTitle)
    % plotPhaseSpace: Plots the velocity vs. x-coordinate phase space
    % for each body in separate subplots and includes a custom figure title.
    %
    % Inputs:
    %   y: Matrix containing the positions and velocities of the bodies.
    %      Each row corresponds to a time step and columns are in the format:
    %      [x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3].
    %   figTitle: A string specifying the title of the figure.

    % Extract the phase space data
    x1 = y(:, 1); vx1 = y(:, 3); % Body 1
    x2 = y(:, 5); vx2 = y(:, 7); % Body 2
    x3 = y(:, 9); vx3 = y(:, 11); % Body 3
    
    % Create the figure
    figure('Name', figTitle, 'NumberTitle', 'off');
    sgtitle(figTitle); % Set the super-title for the figure
    
    % Body 1
    subplot(3, 1, 1);
    plot(x1, vx1, 'r');
    xlabel('x-coordinate'); ylabel('vx');
    title('Phase Space: Body 1');
    grid on;
    
    % Body 2
    subplot(3, 1, 2);
    plot(x2, vx2, 'g');
    xlabel('x-coordinate'); ylabel('vx');
    title('Phase Space: Body 2');
    grid on;
    
    % Body 3
    subplot(3, 1, 3);
    plot(x3, vx3, 'b');
    xlabel('x-coordinate'); ylabel('vx');
    title('Phase Space: Body 3');
    grid on;
end


%Differential Equations
function dydt = threeBodyODE(t, y, masses)
    % Extract positions and velocities from input vector y
    r1 = y(1:2); v1 = y(3:4);
    r2 = y(5:6); v2 = y(7:8);
    r3 = y(9:10); v3 = y(11:12);
    
    % Masses of the three bodies
    m1 = masses(1);
    m2 = masses(2);
    m3 = masses(3);
    
    % Gravitational constant (set to 1 for simplicity; adjust as needed)
    G = 1;
    
    % Compute distances between bodies
    r12 = r2 - r1; d12 = norm(r12);
    r13 = r3 - r1; d13 = norm(r13);
    r23 = r3 - r2; d23 = norm(r23);
    
    % Compute accelerations using Newton's law of gravitation
    a1 = G * ((m2 * r12 / d12^3) + (m3 * r13 / d13^3));
    a2 = G * ((m1 * -r12 / d12^3) + (m3 * r23 / d23^3));
    a3 = G * ((m1 * -r13 / d13^3) + (m2 * -r23 / d23^3));
    
    % Combine derivatives into a single vector
    dydt = [
        v1; % Velocity of body 1
        a1; % Acceleration of body 1
        v2; % Velocity of body 2
        a2; % Acceleration of body 2
        v3; % Velocity of body 3
        a3  % Acceleration of body 3
    ];
end

%Animate function
function animateThreeBody(t, y)
    % Extract position data for each body
    r1 = y(:, 1:2); % Positions of body 1
    r2 = y(:, 5:6); % Positions of body 2
    r3 = y(:, 9:10); % Positions of body 3
    
    % Create the figure
    figure;
    hold on;
    grid on;
    axis equal;
    xlabel('x'); ylabel('y');
    title('Three-Body Problem Animation');
    
    % Define axis limits based on positions
    all_positions = [r1; r2; r3];
    x_min = min(all_positions(:, 1)); x_max = max(all_positions(:, 1));
    y_min = min(all_positions(:, 2)); y_max = max(all_positions(:, 2));
    axis([x_min, x_max, y_min, y_max]);
    
    % Plot initial positions
    body1 = plot(r1(1, 1), r1(1, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    body2 = plot(r2(1, 1), r2(1, 2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    body3 = plot(r3(1, 1), r3(1, 2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    % Plot trajectory traces
    trace1 = plot(r1(1, 1), r1(1, 2), 'r-', 'LineWidth', 1);
    trace2 = plot(r2(1, 1), r2(1, 2), 'g-', 'LineWidth', 1);
    trace3 = plot(r3(1, 1), r3(1, 2), 'b-', 'LineWidth', 1);
    
    % Define the video writer object
    v = VideoWriter('animation.mp4', 'MPEG-4'); % 'MPEG-4' for MP4 format
    v.FrameRate = 30; % Set the frame rate
    open(v); % Open the video file for writing

    % Animation loop
    for i = 1:length(t)
        % Update body positions
        set(body1, 'XData', r1(i, 1), 'YData', r1(i, 2));
        set(body2, 'XData', r2(i, 1), 'YData', r2(i, 2));
        set(body3, 'XData', r3(i, 1), 'YData', r3(i, 2));
        
        % Update trajectory traces
        set(trace1, 'XData', r1(1:i, 1), 'YData', r1(1:i, 2));
        set(trace2, 'XData', r2(1:i, 1), 'YData', r2(1:i, 2));
        set(trace3, 'XData', r3(1:i, 1), 'YData', r3(1:i, 2));
        
        frame = getframe(gcf);
        writeVideo(v, frame); % Write the frame to the video

        % Pause to create animation effect
        pause(0.01);
    end

    % Close the video file
    close(v);

end

%Solvers
function [time, y] = eulerThreeBody(threeBodyODE, masses, y0, tspan, dt)
    % eulerThreeBodyForAnimation: Solves the three-body problem using Euler method
    % and formats output for animation.
    %
    % Inputs:
    %   threeBodyODE - Handle to the ODE function
    %   masses - Masses of the three bodies [m1, m2, m3]
    %   y0 - Initial state vector [r1x, r1y, v1x, v1y, r2x, r2y, v2x, v2y, r3x, r3y, v3x, v3y]
    %   tspan - Time span [t_start, t_end]
    %   dt - Time step
    %
    % Outputs:
    %   time - Time vector
    %   y - Solution matrix for positions and velocities over time
    
    % Time setup
    t_start = tspan(1);
    t_end = tspan(2);
    time = t_start:dt:t_end; % Time vector
    num_steps = length(time);
    
    % Initialize solution matrix
    y = zeros(num_steps, length(y0));
    y(1, :) = y0; % Initial state
    
    % Euler integration loop
    for i = 1:num_steps - 1
        % Current state and time
        t = time(i);
        current_state = y(i, :)';
        
        % Evaluate ODE to get derivatives
        dydt = threeBodyODE(t, current_state, masses);
        
        % Euler update
        y(i + 1, :) = current_state' + dydt' * dt;
    end
end

function [time, y] = rkf45ThreeBody(threeBodyODE, masses, y0, tspan, tol)
    % rkf45ThreeBody: Solves the three-body problem using RKF4(5) adaptive step size.
    %
    % Inputs:
    %   threeBodyODE - Handle to the ODE function
    %   masses - Masses of the three bodies [m1, m2, m3]
    %   y0 - Initial state vector [r1x, r1y, v1x, v1y, r2x, r2y, v2x, v2y, r3x, r3y, v3x, v3y]
    %   tspan - Time span [t_start, t_end]
    %   tol - Tolerance for adaptive step size
    %
    % Outputs:
    %   time - Time vector
    %   y - Solution matrix for positions and velocities over time

    % Coefficients for RKF4(5)
    a = [0, 1/4, 3/8, 12/13, 1, 1/2];
    b = [
        0,         0,         0,       0,         0, 0;
        1/4,       0,         0,       0,         0, 0;
        3/32,      9/32,      0,       0,         0, 0;
        1932/2197, -7200/2197, 7296/2197, 0,      0, 0;
        439/216,   -8,         3680/513, -845/4104, 0, 0;
        -8/27,     2,         -3544/2565, 1859/4104, -11/40, 0;
    ]; 
    c4 = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    c5 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];

    % Initial conditions
    t_start = tspan(1);
    t_end = tspan(2);
    h = (t_end - t_start) / 100; % Initial step size
    t = t_start;
    y = y0(:);
    
    % Store solutions
    time = t;
    solution = y';
    
    % RKF loop
    while t < t_end
        % Ensure the last step ends exactly at t_end
        if t + h > t_end
            h = t_end - t;
        end

        % Compute intermediate stages
        k1 = h * threeBodyODE(t, y, masses);
        k2 = h * threeBodyODE(t + a(2) * h, y + b(2, 1) * k1, masses);
        k3 = h * threeBodyODE(t + a(3) * h, y + b(3, 1) * k1 + b(3, 2) * k2, masses);
        k4 = h * threeBodyODE(t + a(4) * h, y + b(4, 1) * k1 + b(4, 2) * k2 + b(4, 3) * k3, masses);
        k5 = h * threeBodyODE(t + a(5) * h, y + b(5, 1) * k1 + b(5, 2) * k2 + b(5, 3) * k3 + b(5, 4) * k4, masses);
        k6 = h * threeBodyODE(t + a(6) * h, y + b(6, 1) * k1 + b(6, 2) * k2 + b(6, 3) * k3 + b(6, 4) * k4 + b(6, 5) * k5, masses);

        % Compute 4th-order and 5th-order solutions
        y4 = y + c4(1) * k1 + c4(2) * k2 + c4(3) * k3 + c4(4) * k4 + c4(5) * k5 + c4(6) * k6;
        y5 = y + c5(1) * k1 + c5(2) * k2 + c5(3) * k3 + c5(4) * k4 + c5(5) * k5 + c5(6) * k6;

        % Error estimate
        err = norm(y5 - y4, Inf);

        % Accept or reject step based on error
        if err <= tol
            % Accept step
            t = t + h;
            y = y5; % Use the higher-order solution
            time = [time; t];
            solution = [solution; y'];
        end

        % Update step size
        if err == 0
            s = 2; % Avoid division by zero
        else
            s = 0.84 * (tol / err)^(1/4);
        end
        h = min(h * s, t_end - t); % Adjust step size
    end

    % Return solution
    y = solution;
end

function [time, y] = leapfrogThreeBody(threeBodyODE, masses, y0, tspan, dt)
    % leapfrogThreeBody: Solves the three-body problem using the Leapfrog method.
    %
    % Inputs:
    %   threeBodyODE - Handle to the ODE function
    %   masses - Masses of the three bodies [m1, m2, m3]
    %   y0 - Initial state vector [r1x, r1y, v1x, v1y, r2x, r2y, ...]
    %   tspan - Time span [t_start, t_end]
    %   dt - Time step
    %
    % Outputs:
    %   time - Time vector
    %   y - Solution matrix for positions and velocities over time

    % Time setup
    t_start = tspan(1);
    t_end = tspan(2);
    time = t_start:dt:t_end; % Time vector
    num_steps = length(time);
    
    % Initialize state vector
    num_states = length(y0);
    y = zeros(num_steps, num_states); % To store positions and velocities
    y(1, :) = y0; % Initial conditions
    
    % Extract initial positions and velocities
    r1 = y0(1:2); v1 = y0(3:4);
    r2 = y0(5:6); v2 = y0(7:8);
    r3 = y0(9:10); v3 = y0(11:12);
    
    % Initial accelerations
    y_initial = y0(:);
    dydt = threeBodyODE(0, y_initial, masses);
    a1 = dydt(3:4); % Acceleration of body 1
    a2 = dydt(7:8); % Acceleration of body 2
    a3 = dydt(11:12); % Acceleration of body 3
    
    % Leapfrog steps
    for i = 1:num_steps-1
        % Update positions using velocities (half-step)
        r1 = r1 + v1 * dt + 0.5 * a1 * dt^2;
        r2 = r2 + v2 * dt + 0.5 * a2 * dt^2;
        r3 = r3 + v3 * dt + 0.5 * a3 * dt^2;

        % Compute new accelerations
        y_temp = [r1; v1; r2; v2; r3; v3];
        dydt = threeBodyODE(time(i), y_temp, masses);
        a1_new = dydt(3:4); % New acceleration of body 1
        a2_new = dydt(7:8); % New acceleration of body 2
        a3_new = dydt(11:12); % New acceleration of body 3

        % Update velocities using the new accelerations (full-step)
        v1 = v1 + 0.5 * (a1 + a1_new) * dt;
        v2 = v2 + 0.5 * (a2 + a2_new) * dt;
        v3 = v3 + 0.5 * (a3 + a3_new) * dt;

        % Update accelerations for the next step
        a1 = a1_new;
        a2 = a2_new;
        a3 = a3_new;

        % Store the current state
        y(i + 1, :) = [r1; v1; r2; v2; r3; v3];
    end
end

%Truncation Error
function evaluateTruncationError(threeBodyODE, masses, t_ref, y_ref, y0, tspan, dt_values)

    % Interpolate reference solution for comparison
    y_ref_interp = @(t) interp1(t_ref, y_ref, t);
    
    % Evaluate errors for each dt
    global_errors = zeros(size(dt_values));
    for i = 1:length(dt_values)
        dt = dt_values(i);
        [t, y] = eulerThreeBody(threeBodyODE, masses, y0, tspan, dt);
        
        % Compute errors at Euler time points
        errors = vecnorm(y - y_ref_interp(t), 2, 2); % L2 norm of errors
        global_errors(i) = max(errors); % Maximum error (global truncation error)
    end
    
    execution_times = zeros(size(dt_values));
    for i = 1:length(dt_values)
        dt = dt_values(i);
        tic; % Start timer
        eulerThreeBody(@threeBodyODE, masses, y0, tspan, dt);
        execution_times(i) = toc; % Stop timer
    end
    
    % Plot Execution Time vs. Step Size
    figure;
    loglog(dt_values, execution_times, 'o-', 'LineWidth', 2);
    xlabel('Step size (\Delta t)');
    ylabel('Execution Time (s)');
    title('Execution Time vs. Step Size (\Delta t)');
    grid on;


    % Plot truncation errors vs. step sizes
    figure;
    loglog(dt_values, global_errors, 'o-', 'LineWidth', 2);
    xlabel('Step size (\Delta t)');
    ylabel('Global truncation error');
    title('Truncation Error of Euler Method');
    grid on;

    accuracy_per_second = 1 ./ (global_errors .* execution_times);

    % Plot Accuracy per Second vs. Step Size
    figure;
    loglog(dt_values, accuracy_per_second, 'o-', 'LineWidth', 2);
    xlabel('Step size (\Delta t)');
    ylabel('Accuracy per Second');
    title('Accuracy per Second vs. Step Size (\Delta t)');
    grid on;
end

%Energy compute
function energy = computeEnergy(y, masses)
    G = 1; % Gravitational constant
    m1 = masses(1); m2 = masses(2); m3 = masses(3);

    % Extract positions and velocities
    r1 = y(:, 1:2); v1 = y(:, 3:4);
    r2 = y(:, 5:6); v2 = y(:, 7:8);
    r3 = y(:, 9:10); v3 = y(:, 11:12);

    % Kinetic energy
    kinetic = 0.5 * m1 * vecnorm(v1, 2, 2).^2 + ...
              0.5 * m2 * vecnorm(v2, 2, 2).^2 + ...
              0.5 * m3 * vecnorm(v3, 2, 2).^2;

    % Potential energy
    r12 = vecnorm(r1 - r2, 2, 2);
    r13 = vecnorm(r1 - r3, 2, 2);
    r23 = vecnorm(r2 - r3, 2, 2);
    potential = -G * (m1 * m2 ./ r12 + m1 * m3 ./ r13 + m2 * m3 ./ r23);

    % Total energy
    energy = kinetic + potential;
end

function plotPeriodicityComparison(methods, t_values, y_values, y0)
    % plotPeriodicityComparison: Plots periodicity for multiple methods by
    % evaluating the distance to the initial state over time.
    %
    % Inputs:
    %   methods - Cell array of strings with method names (e.g., {'Euler', 'ODE45'})
    %   t_values - Cell array of time vectors for each method
    %   y_values - Cell array of solution matrices for each method
    %   y0 - Initial state vector [x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3]
    %
    % Example:
    %   plotPeriodicityComparison({'Euler', 'ODE45'}, {tEuler, tODE45}, {yEuler, yODE45}, y0)

    % Create the figure
    figure;
    hold on;
    colors = lines(length(methods)); % Distinct colors for each method

    % Loop through each method to calculate and plot periodicity
    for i = 1:length(methods)
        t = t_values{i};
        y = y_values{i};

        % Compute the distance to the initial state
        distances = vecnorm(y - y0', 2, 2); % L2 norm of the difference

        % Plot the periodicity
        plot(t, distances, 'Color', colors(i, :), 'LineWidth', 1.5, ...
            'DisplayName', methods{i});
    end

    % Add labels, title, and legend
    xlabel('Time');
    ylabel('Distance to Initial State');
    title('Periodicity Comparison Across Methods');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
end
