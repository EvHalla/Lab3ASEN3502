%% Task 2.5: Earth-Moon-spacecraft challenge

labDir = fileparts(mfilename('fullpath'));
cd(labDir)

x0 = [0; 13262854; 0; -5482; 5245; 0];
durationHours = 1000;
h = 2;  % RK4 step size [s]

[t, x] = rk4(@three_bp, [0, durationHours*3600], x0, h);
[score, tRef, xRef] = evaluate(t, x, 2, durationHours, ...
    'wody8861@colorado.edu');
fprintf('Task 2.5 three-body challenge score: %.0f hours\n', score)

figure('Color', 'w')
plot(x(:,1)/1e3, x(:,2)/1e3, 'b-', 'LineWidth', 1)
hold on
plot(0, 0, 'ko', 'MarkerFaceColor', 'k')
plot(x(1,1)/1e3, x(1,2)/1e3, 'go', 'MarkerFaceColor', 'g')
plot(x(end,1)/1e3, x(end,2)/1e3, 'ro', 'MarkerFaceColor', 'r')
hold off
axis equal
grid on
xlabel('x position [km]')
ylabel('y position [km]')
title(sprintf('Task 2.5: three-body trajectory, RK4 h = %g s', h))
legend('Spacecraft', 'Earth', 'Start', 'End', 'Location', 'best')
ax = gca;
ax.Toolbar.Visible = 'off';
exportgraphics(gcf, 'lab_2_5_trajectory.png', 'Resolution', 200)

if ~isempty(tRef)
    fprintf('First evaluator failure near %.2f hours; %d reference samples returned.\n', ...
        tRef(1)/3600, size(xRef,1))
end
