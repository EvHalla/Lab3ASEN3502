mu    = 3.986e14;                 
R     = 1e8;                     
vcirc = sqrt(mu / R);             
T     = 2*pi*sqrt(R^3 / mu);      

x0    = [R; 0; 0; 0; vcirc; 0];   
tspan = [0, T];
h     = 100;                      

[tE, xE] = euler(@two_bp, tspan, x0, h);
[tR, xR] = rk4(@two_bp, tspan, x0, h);

figure; hold on; grid on; axis equal;
plot(xE(:,1)/1e3, xE(:,2)/1e3, 'r-', 'LineWidth', 1.2);
plot(xR(:,1)/1e3, xR(:,2)/1e3, 'b-', 'LineWidth', 1.2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');            
plot(R/1e3, 0, 'g^', 'MarkerFaceColor', 'g');        
xlabel('x (km)'); ylabel('y (km)');
title(sprintf('Circular orbit, h = %d s, one period', h));
legend('Euler', 'RK4', 'Earth', 'Start', 'Location', 'best');

closeE = norm(xE(end,1:3)' - x0(1:3));
closeR = norm(xR(end,1:3)' - x0(1:3));
fprintf('Steps taken: %d\n', length(tE) - 1);
fprintf('Euler closure error: %.3e m\n', closeE);
fprintf('RK4   closure error: %.3e m\n', closeR);
fprintf('Function evals: Euler %d, RK4 %d\n', length(tE)-1, 4*(length(tR)-1));