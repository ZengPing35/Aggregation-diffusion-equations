function r = calculate_f
syms x;
syms y;
syms t;
syms xl;
syms yl;
syms C;
syms theta;
syms r;    
%% function
rho(x, y, t) = function_rho_exact(x, y, t);
gad_H = function_vector_grad_H(rho);
% V(x, y) = 0;
V(x, y) = function_V(x, y);
% % W(x, y) = function_W(x, y);
W(x, y) = C;
%% convolution
xl = r * cos(theta);
yl = r * sin(theta);
a=simplify(W(x-xl,y-yl)*rho(xl,yl,t)) * r;
b=int(a,r,0,1);
rhoW=simplify(int(b,theta,0,2*pi));
%% gradient
D = simplify(gad_H + V(x, y) + rhoW);
% D = simplify(gad_H);
Dx = diff(D,x);
Dy = diff(D,y);
c = simplify(rho * Dx);
d = simplify(rho * Dy);
DDx = diff(c, x);
DDy=diff(d, y);
DD = simplify(DDx + DDy);
F = diff(rho, t) - DD;
r = simplify(F);