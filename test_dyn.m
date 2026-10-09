function xdot = test_dyn(~, x)
    a = 1;
    xdot = [x(2); a^2*x(1)];
end