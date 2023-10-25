function r=rho_exact_solution_y_derivative(x,y,t)
% syms x;
% syms y;
% syms t;
% rho(x, y, t) = function_rho_exact(x, y, t);
% r = diff(rho, y);
% r = simplify(r);

r = 2*x^2*y*exp(-t)*(x - 1)^2*(2*y^2 - 3*y + 1);