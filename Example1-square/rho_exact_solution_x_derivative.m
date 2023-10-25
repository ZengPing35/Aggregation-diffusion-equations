function r=rho_exact_solution_x_derivative(x, y, t)
% clear all;
% syms x;
% syms y;
% syms t;
% rho(x, y, t) = function_rho_exact(x, y, t);
% % rho(x, y, t) = exp(-t) * exp(-30 * ((x-0.5)^2 + (y-0.5)^2)) + 1;
% r = diff(rho, x);
% r = simplify(r);


r = 2*x*y^2*exp(-t)*(y - 1)^2*(2*x^2 - 3*x + 1);
