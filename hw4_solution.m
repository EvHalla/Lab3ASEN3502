%% ASEN 3502 Homework 4

clearvars
close all
clc

%% Question 1(a): Explicit Euler update
% $$\mathbf{x}_{i+1}=\mathbf{x}_i+h\mathbf{f}(t_i,\mathbf{x}_i),
% \qquad t_{i+1}=t_i+h.$$

%% Question 1(b): Implementation
% euler is defined below. States are stored as rows.

%% Question 1(c): Harmonic oscillator
[t, x] = euler(@osc_dyn, [0 10], [1; 0], 0.1);

figure('Color','w')
plot(t, x(:,1), 'b-', 'LineWidth', 1.5)
hold on
plot(t, cos(t), 'k--', 'LineWidth', 1.5)
hold off
grid on
xlabel('Time t')
ylabel('First state x_1(t)')
title('Harmonic oscillator: explicit Euler, h = 0.1')
legend('Euler solution', 'Exact: cos(t)', 'Location', 'best')

%% Oscillator interpretation
% Euler's amplitude grows and its phase lags, while cos(t) has unit amplitude.
% At t = 10, the amplitude is 1.645 and the phase lag is 0.0331 rad.
fprintf('Oscillator amplitude at t = 10: %.6f (exact: 1)\n', ...
    norm(x(end,:)))
fprintf('Oscillator phase lag at t = 10: %.6f rad\n', ...
    10 - 100*atan(0.1))

%% Question 2(a): Exact solution
% The characteristic equation r^2 - a^2 = 0 gives roots +/-a:
%
% $$y(t)=C_1e^{at}+C_2e^{-at},\qquad
% y'(t)=aC_1e^{at}-aC_2e^{-at}.$$
%
% The initial conditions give C1 + C2 = 1 and C1 - C2 = 1, so C1 = 1, C2 = 0:
%
% $$y(t)=e^{at}=e^t,\qquad y'(t)=ae^{at}=e^t\quad(a=1).$$

%% Question 2(b): First-order system
% Let x1 = y and x2 = y'. Then
%
% $$\dot{\mathbf{x}}=\left[\begin{array}{c}x_2\\a^2x_1
% \end{array}\right],\qquad
% \mathbf{x}(0)=\left[\begin{array}{c}1\\a\end{array}\right],
% \qquad a=1.$$

%% Question 2(c-d): Verification test
test_ode(@euler, 0.001, 2e-3)
disp('Euler passed the test.')
[t_test, x_test] = euler(@test_dyn, [0 1], [1; 1], 0.001);
fprintf('Maximum absolute error: %.10g; tolerance: %.10g\n', ...
    max(abs(exp(t_test) - x_test(:,1))), 2e-3)

%% Question 3(a): Error trend
N = [10 20 40 80 160 320 640];
h = 1./N;
E_max = zeros(size(h));
for k = 1:numel(h)
    [t_k, x_k] = euler(@test_dyn, [0 1], [1; 1], h(k));
    E_max(k) = max(abs(exp(t_k) - x_k(:,1)));
end

% Fit E_max = C*h^p.
fit_coefficients = polyfit(log(h), log(E_max), 1);
p = fit_coefficients(1);
C = exp(fit_coefficients(2));
fprintf('       N             h             E_max\n')
for k = 1:numel(N)
    fprintf('%8d    %12.8f    %14.8g\n', N(k), h(k), E_max(k))
end
fprintf('Fitted slope p = %.6f; coefficient C = %.6f\n', p, C)

figure('Color','w')
loglog(h, E_max, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
hold on
loglog(h, C*h.^p, 'k--', 'LineWidth', 1.5)
hold off
grid on
xlabel('Step size h')
ylabel('Maximum absolute error E_{max}')
title('Euler convergence for y'''' = y, y(0) = y''(0) = 1')
legend('Measured error', sprintf('Fit: E = %.3f h^{%.3f}', C, p), ...
    'Location', 'best')

%% Question 3(b): Order of convergence
% The slope is 0.982, close to Euler's global order of one.
% Halving h approximately halves the error. At t = 1, x1 = (1+h)^(1/h):
%
% $$E_{\max}(h)=e-(1+h)^{1/h}.$$
%
% Expanding for small h gives
%
% $$E_{\max}(h)=\frac{e}{2}h+O(h^2),$$
close(gcf)

%% Question 3(c): Step size for a maximum error of 1e-6
% The fit gives h = 6.28e-7. The small-h estimate gives h = 2*1e-6/e
% = 7.36e-7, or about 1.36 million steps on [0, 1].
target_error = 1e-6;
h_fit = (target_error/C)^(1/p);
h_asymptotic = 2*target_error/exp(1);
N_estimate = ceil(1/h_asymptotic);
h_estimate = 1/N_estimate;
% Avoid cancellation when calculating the endpoint error.
estimated_error = -exp(1)*expm1(N_estimate*log1p(h_estimate) - 1);
fprintf('Step size from fitted power law: %.8g\n', h_fit)
fprintf('Asymptotic first-order step size: %.8g\n', h_asymptotic)
fprintf('Choose N = %d, h = %.10g; predicted error = %.10g\n', ...
    N_estimate, h_estimate, estimated_error)

%% Local functions
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
