function r=rho_exact_solution_y_derivative(x,y,t)

% R = function_radius;
% r=-(2*y*pi*sin((pi*(x^2 + y^2))/R^2)*exp(-t))/R^2;

r = y*pi*cos((pi*(x^2 + y^2))/2)*exp(-t);

% syms x;
% syms y;
% syms t;
% rho = function_rho_exact(x, y, t);
% r = diff(rho, y)