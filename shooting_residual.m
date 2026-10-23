function F = shooting_residual(rdot0, r0, rf, T, dynamics, h)

% Construct initial six-component state
x0 = [r0(:); rdot0(1); rdot0(2); 0];

% Propagate trajectory using RK4
[t, x] = rk4(dynamics, [0 T], x0, h);

% Final x and y positions
r_final = x(end, 1:2)';

% Target minus actual final position
F = rf(1:2) - r_final;

end