function plot_shooting(x, r0, rf, plotTitle)

figure
plot(0,0,'LineWidth',2)  % TODO: plot trajectory
hold on
plot(0,0,'ko','MarkerFaceColor','k') % TODO: plot initial position
plot(0,0,'rx','MarkerSize',10,'LineWidth',2) % TODO: plot requested final position
plot(0,0,'bo','MarkerSize',8,'LineWidth',2) % TODO: plot actual final position

axis equal
grid on
xlabel('x Position [m]')
ylabel('y Position [m]')
legend('Trajectory','Initial Position','Required Final Position', ...
    'Actual Final Position')
title(plotTitle)

end