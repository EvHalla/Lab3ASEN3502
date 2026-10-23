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
figure;
subplot(2,1,1)
plot(t, errorNorm)
xlabel('Time (hours)'); 
ylabel('||\Delta r|| (m)');
title('Position Difference, RK4 (h = 100 s) vs ode45');
xlim([0,7.2e5]);
subplot(2,1,2)
plot(t,normode45)
xlabel('Time (hours)'); 
ylabel('Geocentric Radius (km)');
title('Geocentric Radius');
xlim([0,7.2e5]);

