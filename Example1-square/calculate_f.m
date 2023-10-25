function r = calculate_f
syms x;
syms y;
syms t;
syms xl;
syms yl;
syms C;
%% function
rho(x, y, t) = function_rho_exact(x, y, t);
gad_H = function_vector_grad_H(rho);
V(x, y) = function_V(x, y);
W(x, y) = function_W(x, y);
% W(x, y) = C;
%% convolution
a = simplify(W(x-xl,y-yl)*rho(xl,yl,t));
b = int(a,xl,0,1);
rhoW = int(b,yl,0,1);
%% gradient
D = simplify(gad_H + V(x, y) + rhoW);
Dx = diff(D,x);
Dy = diff(D,y);
c = simplify(rho * Dx);
d = simplify(rho * Dy);
DDx = diff(c, x);
DDy=diff(d, y);
DD = simplify(DDx + DDy);
F = diff(rho, t) - DD;
r = simplify(F);
