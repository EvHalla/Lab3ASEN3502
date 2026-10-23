
r0 = [7000e3; 0; 0];
rf = [0; 9000e3; 0];
T = 2400;
h = 10;

v0 = [2500; 7200];
f = @(v) shooting_residual(v, r0, rf, T, @two_bp, h);
J = @(v) numjac(f, v, 1);


miss_before = norm(f(v0));

[v_new, info] = newton_sys(f, J, v0, 1, 10);


miss_after = norm(f(v_new));

fprintf('Initial velocity: [%.6f, %.6f] m/s\n', v0);
fprintf('Corrected velocity: [%.6f, %.6f] m/s\n', v_new);
fprintf('Miss before: %.3f km\n', miss_before/1000);
fprintf('Miss after: %.3f km\n', miss_after/1000);

x0 = [r0; v_new(1); v_new(2); 0];
[t, x] = rk4(@two_bp, [0 T], x0, h);

figure;
plot(x(:,1)/1000, x(:,2)/1000, 'b-', 'LineWidth', 1.5);
hold on;
plot(r0(1)/1000, r0(2)/1000, 'go', 'MarkerFaceColor', 'g');
plot(rf(1)/1000, rf(2)/1000, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(x(end,1)/1000, x(end,2)/1000, 'ko', 'MarkerFaceColor', 'k');
axis equal; grid on;
xlabel('x (km)'); ylabel('y (km)');
title('Task 2.4: Iterative Newton Shooting');
legend('Trajectory', 'Start', 'Target', 'Final Position');
saveas(gcf, 'task2_4_converged.png');