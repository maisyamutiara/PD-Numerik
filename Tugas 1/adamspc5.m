function[t, w] = adamspc(f, a, b, n, alpha)
% inisialisasi variabel awal
  h = (b-a)/n;
  t = zeros(n+1, 1);
  w = zeros(n+1, 1);
  t(1) = a;
  w(1) = alpha;

  %mencari nilai awal selanjutya dengan metode runge-kutte
  for i = 1:4
    t(i+1) = t(i) + h;
    m1 = h * f(t(i), w(i));
    m2 = h * f(t(i) + (h/2), w(i) + (m1/2));
    m3 = h * f(t(i) + (h/2), w(i) + (m2/2));
    m4 = h * f(t(i+1), w(i) + m3);
    w(i+1) = w(i) + (m1 + 2*m2 + 2*m3 + m4)/6;
  endfor

  %Algoritma utama
  for i = 5:n
    t(i+1) = t(i) + h;

  % Adams Bashfourth 5-step
    k1 = f(t(i), w(i));
    k2 = f(t(i-1), w(i-1));
    k3 = f(t(i-2), w(i-2));
    k4 = f(t(i-3), w(i-3));
    k5 = f(t(i-4), w(i-4));
    w(i+1) = w(i) + (h/720)*(1901*k1-2774*k2+2616*k3-1274*k4+251*k5);

  % Adams Moulton 4-step
    k0 = f(t(i+1), w(i+1));
    w(i+1) = w(i) + (h/720)*(251*k0 + 646*k1 - 264*k2 + 106*k3 - 19*k4);
  endfor

