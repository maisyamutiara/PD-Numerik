function [t, w] = milspc(f, a, b, n, alpha)
% inisialisasi variabel awal
  h = (b-a)/n;
  t = zeros(n+1, 1);
  w = zeros(n+1, 1);
  t(1) = a;
  w(1) = alpha;

  %mencari nilai awal selanjutya dengan metode runge-kutte
  for i = 1:3
    t(i+1) = t(i) + h;
    k1 = h* f(t(i), w(i));
    k2 = h* f(t(i) + (h/2), w(i) + (k1/2));
    k3 = h* f(t(i) + (h/2), w(i) + (k2/2));
    k4 = h* f(t(i+1), w(i) + k3);
    w(i+1) = w(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
  endfor

  %Algoritma utama
  for i = 4:n
    t(i+1) = t(i) + h;

  %Predictor Metode Milne
    m1 = f(t(i), w(i));
    m2 = f(t(i-1), w(i-1));
    m3 = f(t(i-2), w(i-2));
    m4 = f(t(i-3), w(i-3));
    w(i+1) = w(i-3) + (4*h/3) * (2*m1 - m2 + 2*m3);

  %Corrector Metode Simpson
    m0 = f(t(i+1), w(i+1));
    w(i+1) = w(i-1) + (h/3)*(m0 + 4*m1 + m2);
  endfor
