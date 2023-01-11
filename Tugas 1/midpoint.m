function [t, w] = midpoint(f, a, b, n, alpha)
  h = (b-a)/n;
  t = zeros(n+1, 1);
  w = zeros(n+1, 1);
  t(1) = a;
  w(1) = alpha;
  for i = 1:n
    t(i+1) = t(i) + h;
    m1 = f(t(i), w(i));
    m2 = f(t(i) + (h/2), w(i)+ (h/2) *m1);
    w(i+1) = w(i) + h*m2;
  endfor
endfunction
