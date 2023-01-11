f = @(t, y) (1+(t-y)^2);
a = 2;
b = 3;
n = 5;
alpha = 1;
[t1, w1] = milspc(f, a, b, n, alpha);

wex = [];
sln = @(t,y) (t + (1-t)^(-1));
for i = 1:length(t1)
  wex(i) = sln(t1(i));
endfor

format long
[t1, w1, wex']

fplot(sln, [2, 3], 'b');
hold on;
scatter(t1, w1, 'r');
legend('Eksak', 'Milne-Simpson');
title('Eksak vs Milne-Simpson Predictor Corrector');
