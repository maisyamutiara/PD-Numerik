f = @(t, y) (y-(t^2)+1);
a = 0;
b = 2;
n = 20;
alpha = 0.5;
[t1, w1] = adamspc5(f, a, b, n, alpha);
[t2, w2] = adamspc(f, a, b, n, alpha);

wex = [];
sln = @(t,y) ((t+1)^2-(0.5*exp(t)));
for i = 1:length(t1)
  wex(i) = sln(t1(i));
endfor

format long
[t1, w1, w2, wex']

fplot(sln, [0, 2], 'b');
hold on;
scatter(t1, w1, 'r');
scatter(t2, w2, 'g');
legend('Eksak', 'Adams PC 5 step', 'Adams PC 4 step');
title('Adams Predictor-Corrector 4 Step vs 5 Step');
