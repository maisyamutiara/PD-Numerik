function[t, w] = euler(f, a, b, n, alpha)
  h = (b - a)/n;
  t = zeros(n+1, 1);
  w = zeros(n+1, 1);
  t(1) = a;
  w(1) = alpha;
  for i = 1: n
    t(i+1) = t(i) + h;
    m1 = f(t(i), w(i));
    w(i+1) = w(i) + h *m1;
  endfor
endfunction