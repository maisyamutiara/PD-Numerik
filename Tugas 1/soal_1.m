f = @(t, y) (y^2/(1+t));
a = 1;
b = 2;
n = 20;
alpha = -log(2)^(-1);
[t1, w1] = euler(f, a, b, n, alpha);
[t2, w2] = midpoint(f, a, b, n, alpha);
[t3, w3] = rk4(f, a, b, n, alpha);

wex = [];
sln = @(t,y) (-(log(t+1))^(-1));
for i = 1:length(t2)
  wex(i) = sln(t2(i));
endfor

format long
[w1, w2, w3, wex', t1]

fplot(sln, [1, 2], 'b');
hold on;
scatter(t1, w1, 'r');
scatter(t2, w2, 'g');
scatter(t3, w3, 'c');
legend('Fungsi Eksak', 'Euler', 'Midpoint', 'Runge-Kutta orde 4');
legend("location", "northwest");
title('Eksak vs Euler vs Midpoint vs Runge-Kutta 4');
