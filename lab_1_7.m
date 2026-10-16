clear all, close all;


h = 10;
x0 = [0; 13262854; 0; -5482; 5245; 0];
tspan = [0, 200*60^2];

[t, x] = rk4(@two_bp, tspan, x0, h);

[score, tRef, xRef] = evaluate(t, x, 1, 200, 'wody8861@colorado.edu');