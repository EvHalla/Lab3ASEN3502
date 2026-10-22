%initial velocity guess for shooting

clear; 
clc; 
close all;

r0 = [7000e3; 0; 0];  
rf = [0; 9000e3; 0];  
T = 2400;  h = 10;
v_guess = [2500; 7200];  %m/s
[~, x] = rk4(@two_bp, [0 T], [r0; v_guess; 0], h);

fprintf('Guess [%g %g] m/s, miss = %.1f km\n', v_guess, norm(x(end,1:3).' - rf)/1e3);
plot_shooting(x, r0, rf, 'Initial Guess');
