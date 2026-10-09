function [t, x] = euler(dynamics, tspan, x0, h)
    validateattributes(h, {'numeric'}, {'scalar','real','finite','positive'})
    validateattributes(tspan, {'numeric'}, {'vector','numel',2,'real','finite'})
    validateattributes(x0, {'numeric'}, {'vector','nonempty','finite'})
    assert(tspan(2) >= tspan(1), 'euler:InvalidInterval', ...
        'The final time must be greater than or equal to the initial time.')
    t = (tspan(1):h:tspan(2)).';
    if t(end) < tspan(2)
        t(end+1,1) = tspan(2);
    end
    x = zeros(numel(t), numel(x0));
    x(1,:) = x0(:).';
    for i = 1:numel(t)-1
        dt = t(i+1) - t(i);
        xdot = dynamics(t(i), x(i,:).');
        assert(iscolumn(xdot) && numel(xdot) == numel(x0), ...
            'euler:DerivativeSize', 'dynamics must return a state-sized column vector.')
        x(i+1,:) = x(i,:) + dt*xdot.';
    end
end

function xdot = osc_dyn(~, x)
    xdot = [x(2); -x(1)];
end

function xdot = test_dyn(~, x)
    a = 1;
    xdot = [x(2); a^2*x(1)];
end

function test_ode(integrator, h, atol)
    [t, x] = integrator(@test_dyn, [0 1], [1; 1], h);
    max_error = max(abs(exp(t) - x(:,1)));
    assert(max_error < atol, '%s failed: max error %g is not below %g', ...
        func2str(integrator), max_error, atol)
end
