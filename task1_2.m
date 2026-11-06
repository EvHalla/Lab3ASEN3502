%% Task 1.2: Euler circular-orbit trajectory and specific orbital energy
% Uses the existing dynamics and integrator without modifying either.
labDir = fileparts(mfilename('fullpath'));
mu = 3.986e14;                 % m^3/s^2
R = 1e8;                      % m
vcirc = sqrt(mu/R);            % m/s
T = 2*pi*sqrt(R^3/mu);         % s
h = 10;                       % s
x0 = [R; 0; 0; 0; vcirc; 0];
[t, x] = euler(@two_bp, [0 T], x0, h);

radius = sqrt(sum(x(:,1:3).^2, 2));
speedSquared = sum(x(:,4:6).^2, 2);
energy = 0.5*speedSquared - mu./radius; % J/kg
initialEnergy = -mu/(2*R);

fig = figure('Color', 'w', 'Theme', 'light', 'Position', [100 100 1200 520]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
theta = linspace(0, 2*pi, 1000);
plot(R*cos(theta)/1e3, R*sin(theta)/1e3, 'k--', 'LineWidth', 1)
hold on
plot(x(:,1)/1e3, x(:,2)/1e3, 'b-', 'LineWidth', 1.6)
plot(0, 0, 'ko', 'MarkerFaceColor', 'k')
plot(x(1,1)/1e3, x(1,2)/1e3, 'g^', 'MarkerFaceColor', 'g')
plot(x(end,1)/1e3, x(end,2)/1e3, 'ro', 'MarkerFaceColor', 'r')
hold off
axis equal
grid on
xlabel('x position (km)')
ylabel('y position (km)')
title('Trajectory over one circular-orbit period')
legend('Exact circular orbit', 'Euler, h = 10 s', 'Earth', 'Start', 'End', ...
    'Location', 'southoutside', 'NumColumns', 2)

nexttile
plot(t/3600, energy/1e6, 'b-', 'LineWidth', 1.6)
hold on
yline(initialEnergy/1e6, 'k--', 'LineWidth', 1)
hold off
grid on
xlim([0 T/3600])
energyRange = (max(energy)-min(energy))/1e6;
ylim([min(energy)/1e6 - 0.08*energyRange, max(energy)/1e6 + 0.08*energyRange])
xlabel('Time (hours)')
ylabel('Specific orbital energy (MJ/kg)')
title('Specific orbital energy over time')
legend('Euler, h = 10 s', 'Exact constant energy', 'Location', 'southoutside')
sgtitle(sprintf('Task 1.2: R = 10^8 m, Euler h = 10 s, T = %.2f hours', T/3600))

exportgraphics(fig, fullfile(labDir, 'task1_2.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white')
fprintf('Task 1.2: T = %.8f s; speed = %.8f m/s\n', T, vcirc)
fprintf('Initial energy = %.9g J/kg; final energy = %.9g J/kg\n', ...
    energy(1), energy(end))
fprintf('Energy increase = %.9g J/kg (%.6f%% of initial magnitude)\n', ...
    energy(end)-energy(1), 100*(energy(end)-energy(1))/abs(energy(1)))
