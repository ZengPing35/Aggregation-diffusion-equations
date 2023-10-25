function rho = function_rho_exact(x, y, t)

rho = exp(-t) * x^2* y^2 * (1-x)^2 * (1-y)^2 + 0.1;
