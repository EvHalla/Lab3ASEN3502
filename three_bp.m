function xdot = three_bp(t, x)

    muE = 3.986e14;       % Earth's gravitational parameter [m^3/s^2]
    muM = 4.904e12;       % Moon's gravitational parameter [m^3/s^2]
    RM = 3.844e8;        % Prescribed Moon orbital radius [m]
    theta0 = 2.3243;      % Moon's initial angle [rad]

    r = x(1:3);
    v = x(4:6);

    nM = sqrt(muE / RM^3);
    thetaM = nM*t + theta0;
    rM = RM * [cos(thetaM); sin(thetaM); 0];
    dM = rM - r;

    aE = -muE * r / norm(r)^3;
    aM = muM * dM / norm(dM)^3;
    xdot = [v; aE + aM];
end
