clc;
clear;
close all
x0new = [0; 13262854; 0; -5482; 5245; 0];
intnew = [0, 200*3600];
h=100;
[t, x] = rk4(@two_bp, intnew, x0new, h);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[t2,x2] = ode45(@two_bp, intnew, x0new,options);
x_ode45 = interp1(t2, x2, t);
for i = 1:length(x)
errorNorm(i,1) = norm(x(i,1:3) - x_ode45(i,1:3));
normode45(i,1) = norm(x_ode45(i,1:3));
end

%overlay of both trajectories
figure; 
hold on; 
grid on;
plot(x(:,1)/1e3,  x(:,2)/1e3,  'b-',  'LineWidth', 1.2);
plot(x_ode45(:,1)/1e3, x_ode45(:,2)/1e3, '--g', 'LineWidth', 1.2);
xlabel('x (km)'); ylabel('y (km)');
legend('RK4, h = 100 s', 'ode45, RelTol 1e-12');
title('Trajectory Overlay');

%position difference and geocentric radius versus time plots
fig5 = figure('Color', 'w', 'Theme', 'light', 'Position', [100 100 1100 750]);
tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
plot(t/3600, errorNorm, 'b-', 'LineWidth', 1.4)
grid on
xlabel('Time (hours)'); 
ylabel('||\Delta r|| (m)');
title('Position Difference, RK4 (h = 100 s) vs ode45');
xlim([0,intnew(2)/3600]);
nexttile
plot(t/3600,normode45/1e3, 'b-', 'LineWidth', 1.4)
grid on
xlabel('Time (hours)'); 
ylabel('Geocentric Radius (km)');
title('Geocentric Radius from ode45');
xlim([0,intnew(2)/3600]);
sgtitle({'Task 2.1: fixed-step RK4 versus adaptive ode45', ...
    'RK4 h = 100 s; ode45 RelTol = 10^{-12}, AbsTol = 10^{-14}'})
labDir = fileparts(mfilename('fullpath'));
exportgraphics(fig5, fullfile(labDir, 'task2_1_figure5.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white')
[peakError, peakIndex] = max(errorNorm);
fprintf('Peak position difference: %.6f m at %.6f hours\n', ...
    peakError, t(peakIndex)/3600)
fprintf('Minimum sampled geocentric radius: %.6f km\n', min(normode45)/1e3)

