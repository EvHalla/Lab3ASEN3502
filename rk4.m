function [t, x] = rk4(dynamics, tspan, x0, h)

    % Initialize time vector
    t = (tspan(1):h:tspan(2))';

    % Include final time if necessary
    if t(end) < tspan(2)
        t(end+1) = tspan(2);
    end

    % Initialize state matrix
    x0 = x0(:);
    x = zeros(length(t), length(x0));
    x(1,:) = x0';

    % RK4 integration
    for n = 1:length(t)-1

        tn = t(n);
        xn = x(n,:)';

        dt = t(n+1) - t(n);

        k1 = dynamics(tn, xn);
        k2 = dynamics(tn + dt/2, xn + dt*k1/2);
        k3 = dynamics(tn + dt/2, xn + dt*k2/2);
        k4 = dynamics(tn + dt, xn + dt*k3);

        x(n+1,:) = (xn + dt/6*(k1 + 2*k2 + 2*k3 + k4))';
    end
end