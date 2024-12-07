clear; clc; close all;
G = 4 * pi^2;
% Initial conditions
masses = [1, 0.0009546*30, 0.0002858]; % Sun, Jupiter, Saturn
% Jupiter's orbital parameters
a_jupiter = 5.2; % Semi-major axis (AU) Default: 5.2
e_jupiter = 0.0489; % Eccentricity Default: 0.0489
r_jupiter = a_jupiter * (1 - e_jupiter); % Perihelion distance
v_jupiter = sqrt(G * masses(1) * (2 / r_jupiter - 1 / a_jupiter)); % Velocity at perihelion

% Saturn's orbital parameters
a_saturn = 9.5; % Semi-major axis (AU) Default: 9.5
e_saturn = 0.0520; % Eccentricity Default = 0.0520
r_saturn = a_saturn * (1 - e_saturn); % Perihelion distance
v_saturn = sqrt(G * masses(1) * (2 / r_saturn - 1 / a_saturn)); % Velocity at perihelion

% Initial conditions
y0 = [
    0;         0;           % Initial position of the Sun
    0;         0;           % Initial velocity of the Sun
    r_jupiter; 0;           % Initial position of Jupiter (perihelion)
    0;         v_jupiter;   % Initial velocity of Jupiter (tangential at perihelion)
    r_saturn;  0;           % Initial position of Saturn (perihelion)
    0;         v_saturn     % Initial velocity of Saturn (tangential at perihelion)
];
tspan = [0 30000]; % Time range for the simulation 50000 for mass experiments
dt = 0.01;
dtLF = 0.01;

plotInitialConditions(y0);

% Solve using Leapfrog
tic;
[timeLF, yLF] = leapfrogThreeBody(@threeBodyODE, masses, y0, tspan, dtLF);
computation_times(4) = toc;
plotPhaseSpace(yLF, "Leapfrog");

r1LF = yLF(:, 1:2);
r2LF = yLF(:, 5:6);
r3LF = yLF(:, 9:10);
figure;
a = 320000;
b = 350000;
plot(r1LF(a:b, 1), r1LF(a:b, 2), 'ro', r2LF(a:b, 1), r2LF(a:b, 2), 'g', r3LF(a:b, 1), r3LF(a:b, 2), 'b');
xlabel('x'); ylabel('y'); legend('Sun', 'Jupiter', 'Saturn');
title('Partial Trajectory');

% Visualise the Leapfrog method solution (e.g., positions of the three bodies)
r1LF = yLF(:, 1:2);
r2LF = yLF(:, 5:6);
r3LF = yLF(:, 9:10);
figure;
plot(r1LF(:, 1), r1LF(:, 2), 'ro', r2LF(:, 1), r2LF(:, 2), 'g', r3LF(:, 1), r3LF(:, 2), 'b');
xlabel('x'); ylabel('y'); legend('Sun', 'Jupiter', 'Saturn');
title('Sun-Jupiter-Saturn Trajectories (Leapfrog)');

%Truncation error evaluation (Euler)
dt_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005];

periodic_differences = zeros(size(dt_values));
%{
threshold = 1000; % Set a radial distance threshold (e.g., 100 AU)
unboundedTime = detectUnboundedTime(timeLF, yLF, masses, threshold);
%}
% Plot periodic differences
figure;
plot(dt_values, periodic_differences, '-o', 'LineWidth', 2);
xlabel('Step size (\Delta t)');
ylabel('Periodic Difference (\|y(T) - y(0)\|)');
title('Periodicity Check for Different \Delta t');
grid on;

energyLF = computeEnergy(yLF, masses);
% Plot energy conservation for leapfrog
figure;
plot(timeLF, energyLF, 'g-.', 'LineWidth', 1.5); hold on;
xlabel('Time');
ylabel('Total Energy');
title('Energy Conservation Comparison for Leapfrog');
legend('Leapfrog');
grid on;

periods = calculateOrbitalPeriods(timeLF, yLF, masses);

unboundedTime = detectUnboundedEnergy(timeLF, yLF, masses, G);

v = VideoWriter('animation.mp4', 'MPEG-4'); % 'MPEG-4' for MP4 format
v.FrameRate = 30; % Set the frame rate
open(v); % Open the video file for writing
animateThreeBody(timeLF, yLF);
% Close the video file
close(v);

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


% Differential Equations
function dydt = threeBodyODE(t, y, masses)
    % Extract positions and velocities from input vector y
    r2 = y(5:6); v2 = y(7:8); % Position and velocity of body 2 (Jupiter)
    r3 = y(9:10); v3 = y(11:12); % Position and velocity of body 3 (Saturn)
    
    % Masses of the bodies
    m1 = masses(1); % Sun
    m2 = masses(2); % Jupiter
    m3 = masses(3); % Saturn
    
    % Gravitational constant (set to 1 for simplicity; adjust as needed)
    G = 4 * pi^2;
    
    % Compute distances between bodies
    r12 = r2; d12 = norm(r12); % Distance between Sun and Jupiter
    r13 = r3; d13 = norm(r13); % Distance between Sun and Saturn
    r23 = r3 - r2; d23 = norm(r23); % Distance between Jupiter and Saturn
    
    % Compute accelerations
    % The Sun (Body 1) is stationary and exerts forces but does not move
    a2 = G * (-m1 * r12 / d12^3 + m3 * r23 / d23^3); % Jupiter
    a3 = G * (-m1 * r13 / d13^3 - m2 * -r23 / d23^3); % Saturn
    
    % Combine derivatives into a single vector
    dydt = [
        0; 0;               % Velocity of the Sun (fixed)
        0; 0;               % Acceleration of the Sun (fixed)
        v2; a2;             % Velocity and acceleration of Jupiter
        v3; a3              % Velocity and acceleration of Saturn
    ];
end


%Animate function
%{
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
    title('Sun-Jupiter-Saturn Animation');
    
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
        
        % Pause to create animation effect
        pause(0.01);
    end
end
%}

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
    title('Sun-Jupiter-Saturn Animation');
    
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
    
    % Frame skipping to speed up animation
    skip_frames = 1000; % Plot every 10th frame
    
    % Disable axis updates for better performance
    set(gcf, 'GraphicsSmoothing', 'off');
    
    % Define the video writer object
    v = VideoWriter('animation.mp4', 'MPEG-4'); % 'MPEG-4' for MP4 format
    v.FrameRate = 30; % Set the frame rate
    open(v); % Open the video file for writing

    % Animation loop
    for i = 1:skip_frames:length(t)
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

        % Reduce pause duration for faster updates
        pause(0.0001);
    end

    % Close the video file
    close(v);
end


%Solvers
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
    G = 4 * pi^2; % Gravitational constant
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

function periods = calculateOrbitalPeriods(time, y, masses)
    % calculateOrbitalPeriods: Calculates the orbital period of planets
    % based on their radial distance from the Sun.
    %
    % Inputs:
    %   time - Time vector from the simulation
    %   y - Solution matrix from the simulation
    %   masses - Masses of the bodies (for labeling purposes)
    %
    % Outputs:
    %   periods - A struct containing the orbital periods of the planets

    % Extract radial distances of Jupiter and Saturn from the Sun
    r2 = y(:, 5:6); % Jupiter's position
    r3 = y(:, 9:10); % Saturn's position
    r2_dist = vecnorm(r2, 2, 2); % Radial distance of Jupiter from Sun
    r3_dist = vecnorm(r3, 2, 2); % Radial distance of Saturn from Sun

    % Find peaks in the radial distance to calculate orbital periods
    [~, locs2] = findpeaks(r2_dist, time); % Peaks for Jupiter
    [~, locs3] = findpeaks(r3_dist, time); % Peaks for Saturn

    % Calculate periods as the difference between consecutive peaks
    period_jupiter = mean(diff(locs2)); % Average time between peaks
    period_saturn = mean(diff(locs3)); % Average time between peaks

    % Store results in a struct
    periods = struct('Jupiter', period_jupiter, 'Saturn', period_saturn);

    % Display results
    fprintf('Orbital Period of Jupiter: %.3f years\n', period_jupiter);
    fprintf('Orbital Period of Saturn: %.3f years\n', period_saturn);

    dt = mean(diff(time)); % Time step
    Fs = 1 / dt;           % Sampling frequency

    % Compute the FFT of the radial distances
    n = length(r2_dist);       % Number of data points
    n1 = length(r3_dist);
    fft_result = fft(r2_dist); % Perform FFT
    fft_result1 = fft(r3_dist);
    fft_magnitude = abs(fft_result / n); % Normalized magnitude
    fft_magnitude1 = abs(fft_result1 / n1); % Normalized magnitude
    fft_magnitude = fft_magnitude(1:n/2+1); % One-sided spectrum
    fft_magnitude1 = fft_magnitude1(1:n1/2+1); % One-sided spectrum
    fft_magnitude(2:end-1) = 2 * fft_magnitude(2:end-1); % Double-sided to one-sided scaling
    fft_magnitude1(2:end-1) = 2 * fft_magnitude1(2:end-1); % Double-sided to one-sided scaling

    % Frequency vector for the FFT
    f = Fs * (0:(n/2)) / n;
    f1 = Fs * (0:(n1/2)) / n1;

    % Plot the FFT results
    figure;
    hold on;
    plot(f, fft_magnitude, 'LineWidth', 1.5, 'Color', 'g');
    plot(f1,fft_magnitude1, 'LineWidth',1.5, 'Color', 'b')
    xlabel('Frequency (1/Time Unit)');
    ylabel('Magnitude');
    title('FFT of Radial Distances');
    grid on;
    hold off;

    % Plot radial distances and detected peaks
    figure;
    hold on;
    plot(time, r2_dist, 'g', 'LineWidth', 1.5, 'DisplayName', 'Jupiter Distance');
    plot(time, r3_dist, 'b', 'LineWidth', 1.5, 'DisplayName', 'Saturn Distance');
    xlabel('Time');
    ylabel('Radial Distance (AU)');
    title('Radial Distance from Sun and Orbital Periods');
    legend('show', 'Location', 'best');
    grid on;
    hold off;



end

%{
function unboundedTime = detectUnboundedTime(time, y, masses, threshold)
    % detectUnboundedTime: Determines the time it takes for a body to become unbounded.
    %
    % Inputs:
    %   time - Time vector from the simulation
    %   y - Solution matrix from the simulation (positions and velocities)
    %   masses - Masses of the Sun, Jupiter, and Saturn (for labeling purposes)
    %   threshold - Radial distance threshold (e.g., 100 AU) to consider unbounded
    %
    % Outputs:
    %   unboundedTime - Struct containing the time each body becomes unbounded
    %                   If a body does not become unbounded, the value is NaN.

    % Extract radial distances of Jupiter and Saturn from the Sun
    r2 = vecnorm(y(:, 5:6), 2, 2); % Jupiter's radial distance
    r3 = vecnorm(y(:, 9:10), 2, 2); % Saturn's radial distance

    % Find the first time step where each body exceeds the threshold
    unboundedTime.Jupiter = findUnboundedTime(time, r2, threshold);
    unboundedTime.Saturn = findUnboundedTime(time, r3, threshold);

    % Display results
    if isnan(unboundedTime.Jupiter)
        fprintf('Jupiter did not become unbounded within the simulation time.\n');
    else
        fprintf('Jupiter became unbounded at t = %.3f years.\n', unboundedTime.Jupiter);
    end

    if isnan(unboundedTime.Saturn)
        fprintf('Saturn did not become unbounded within the simulation time.\n');
    else
        fprintf('Saturn became unbounded at t = %.3f years.\n', unboundedTime.Saturn);
    end
end

function t_unbounded = findUnboundedTime(time, r, threshold)
    % findUnboundedTime: Finds the first time a radial distance exceeds a threshold.
    %
    % Inputs:
    %   time - Time vector from the simulation
    %   r - Radial distance vector
    %   threshold - Radial distance threshold
    %
    % Outputs:
    %   t_unbounded - Time when the radial distance first exceeds the threshold

    % Find the index of the first instance where the radial distance exceeds the threshold
    idx = find(r > threshold, 1, 'first');

    % If no index is found, return NaN
    if isempty(idx)
        t_unbounded = NaN;
    else
        t_unbounded = time(idx);
    end
end
%}

function unboundedTime = detectUnboundedEnergy(time, y, masses, G)
    % detectUnboundedEnergy: Determines the time it takes for a body to become unbounded
    % based on its total energy (kinetic + potential).
    %
    % Inputs:
    %   time - Time vector from the simulation
    %   y - Solution matrix from the simulation (positions and velocities)
    %   masses - Masses of the Sun, Jupiter, and Saturn [m1, m2, m3]
    %   G - Gravitational constant
    %
    % Outputs:
    %   unboundedTime - Struct containing the time each body becomes unbounded
    %                   If a body does not become unbounded, the value is NaN.

    % Extract positions and velocities
    r2 = y(:, 5:6); % Jupiter's position
    r3 = y(:, 9:10); % Saturn's position
    v2 = y(:, 7:8); % Jupiter's velocity
    v3 = y(:, 11:12); % Saturn's velocity

    % Masses
    m1 = masses(1); % Sun
    m2 = masses(2); % Jupiter
    m3 = masses(3); % Saturn

    % Calculate total energy for each body at each time step
    energyJupiter = computeTotalEnergy(r2, v2, m1, m2, G);
    energySaturn = computeTotalEnergy(r3, v3, m1, m3, G);

    % Smooth energy data to prevent false positives
    energyJupiter = smoothdata(energyJupiter, 'movmean', 10);
    energySaturn = smoothdata(energySaturn, 'movmean', 10);

    % Find the first time step where total energy becomes positive
    unboundedTime.Jupiter = findUnboundedEnergy(time, energyJupiter);
    unboundedTime.Saturn = findUnboundedEnergy(time, energySaturn);

    % Display results
    if isnan(unboundedTime.Jupiter)
        fprintf('Jupiter did not become unbounded within the simulation time.\n');
    else
        fprintf('Jupiter became unbounded at t = %.3f years.\n', unboundedTime.Jupiter);
    end

    if isnan(unboundedTime.Saturn)
        fprintf('Saturn did not become unbounded within the simulation time.\n');
    else
        fprintf('Saturn became unbounded at t = %.3f years.\n', unboundedTime.Saturn);
    end
end

function totalEnergy = computeTotalEnergy(r, v, mSun, mPlanet, G)
    % computeTotalEnergy: Computes the total energy (kinetic + potential) of a planet.
    %
    % Inputs:
    %   r - Position vector of the planet (Nx2 array)
    %   v - Velocity vector of the planet (Nx2 array)
    %   mSun - Mass of the Sun
    %   mPlanet - Mass of the planet
    %   G - Gravitational constant
    %
    % Outputs:
    %   totalEnergy - Total energy of the planet at each time step (Nx1 array)

    % Radial distance from the Sun
    r_dist = vecnorm(r, 2, 2);

    % Kinetic energy: (1/2) * m * v^2
    kineticEnergy = 0.5 * mPlanet * vecnorm(v, 2, 2).^2;

    % Potential energy: -G * m1 * m2 / r
    potentialEnergy = -G * mSun * mPlanet ./ r_dist;

    % Total energy
    totalEnergy = kineticEnergy + potentialEnergy;
end

function t_unbounded = findUnboundedEnergy(time, energy)
    % findUnboundedEnergy: Finds the first time total energy becomes positive.
    %
    % Inputs:
    %   time - Time vector from the simulation
    %   energy - Total energy of the body at each time step
    %
    % Outputs:
    %   t_unbounded - Time when total energy first becomes positive

    % Ignore initial time steps to avoid numerical artifacts
    initialBuffer = 10; % Ignore the first 10 time steps
    energy = energy(initialBuffer:end);
    time = time(initialBuffer:end);

    % Find the index where total energy first becomes positive
    idx = find(energy > 0, 1, 'first');

    % If no index is found, return NaN
    if isempty(idx)
        t_unbounded = NaN;
    else
        t_unbounded = time(idx);
    end
end
