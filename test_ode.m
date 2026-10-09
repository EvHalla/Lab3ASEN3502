function test_ode(integrator, h, atol)
    [t, x] = integrator(@test_dyn, [0 1], [1; 1], h);
    max_error = max(abs(exp(t) - x(:,1)));
    assert(max_error < atol, '%s failed: max error %g is not below %g', ...
        func2str(integrator), max_error, atol)
end