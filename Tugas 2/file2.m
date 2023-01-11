f = @(x, y, yp) ((1/8)*(32 + 2*x^3 - y*yp));
fy = @(x, y, yp) (-yp/8);
fyp = @(x, y, yp) (-y/8);
a = 1;
b = 3;
n = 20;
alpha = 17;
beta = 43/3;
m = 10;
tol = 10^(-5);

[xi, w1i, w2i] = nonlinshoot(f, fy, fyp, a, b, n, alpha, beta, 20, tol);

sln = @(x) ((x^2) + (16/x));
w = [];
for i = 1:length(xi)
  w(i) = sln(xi(i));
endfor

[xi', w1i', w']

hold on;
fplot(sln, [1,3], 'k');
scatter(xi, w1i, 'r');
legend('Eksak', 'Aproksimasi');
