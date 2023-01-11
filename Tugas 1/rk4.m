function[t, w] = rk4(f, a, b, n, alpha)
  h = (b - a)/n;
  t = zeros(n+1, 1);
  w = zeros(n+1, 1);
  t(1) = a;
  w(1) = alpha;
  for i = 1: n
    t(i+1) = t(i) + h;
    k1 = h * f(t(i), w(i));
    k2 = h * f(t(i) + (h/2), w(i) + (k1/2));
    k3 = h * f(t(i) + (h/2), w(i) + (k2/2));
    k4 = h * f(t(i+1), w(i) + k3); 
    w(i+1) = w(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
  endfor
endfunction