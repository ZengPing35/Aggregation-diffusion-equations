function r = circumcenter(x1, y1, x2, y2, x3, y3)

d1 = (x2^2 + y2^2) - (x1^2 + y1^2);
d2 = (x3^2 + y3^2) - (x2^2 + y2^2);
d = 2 * ((y3 - y2) * (x2 - x1) - (y2 - y1) * (x3 - x2));
x = ((y3 - y2) * d1 - (y2 - y1) * d2) / d;
y = ((x2 - x1) * d2 - (x3 - x2) * d1) / d;
r = [x; y];
