function result=function_f_rho(x, y, t)

result = -(exp(-t)*(100*sin((pi*(x^2 + y^2))/2) + 220*pi*cos((pi*(x^2 + y^2))/2) + pi*sin(pi*(x^2 + y^2)) + x^2*pi^2*cos(pi*(x^2 + y^2)) + y^2*pi^2*cos(pi*(x^2 + y^2)) - 110*x^2*pi^2*sin((pi*(x^2 + y^2))/2) - 110*y^2*pi^2*sin((pi*(x^2 + y^2))/2) + 1000))/100;
