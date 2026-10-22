function plot_shooting(x, r0, rf, plotTitle)

figure
plot(x(:,1), x(:,2), 'LineWidth',2)  % TODO: plot trajectory
hold on
plot(r0(1),r0(2),'ko','MarkerFaceColor','k') % TODO: plot initial position
plot(rf(1),rf(2),'rx','MarkerSize',10,'LineWidth',2) % TODO: plot requested final position
plot(x(end,1),x(end,2),'bo','MarkerSize',8,'LineWidth',2) % TODO: plot actual final position

axis equal
grid on
xlabel('x Position [m]')
ylabel('y Position [m]')
legend('Trajectory','Initial Position','Required Final Position', ...
    'Actual Final Position')
title(plotTitle)

end
