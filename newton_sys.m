function [x, info] = newton_sys(f, J, x0, atol, maxit)
x = x0;
info = struct('iterations', 0, 'converged', false, 'residualNorm', norm(f(x)));

for iteration = 1:maxit
    residual = f(x);
    if norm(residual) <= atol
        info.iterations = iteration - 1;
        info.converged = true;
        info.residualNorm = norm(residual);
        return
    end
    dx = gauss_pivot(J(x), -residual);

    x = x + dx;
end

info.iterations = maxit;
info.residualNorm = norm(f(x));
info.converged = info.residualNorm <= atol;
