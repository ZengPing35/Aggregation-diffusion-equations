function rho = function_rho_exact(x, y, t)

rho = exp(-t) * sin(0.5 * pi * (x^2+ y^2)) + 10*exp(-t);
