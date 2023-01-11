function[t, w] = adamspc(f, a, b, n, alpha)
% inisialisasi variabel awal
  h = (b-a)/n;
  t = zeros(n+1, 1);
  w = zeros(n+1, 1);
  t(1) = a;
  w(1) = alpha;

  %mencari nilai awal selanjutya dengan metode runge-kutte
  for i = 1:3
    t(i+1) = t(i) + h;
    m1 = h * f(t(i), w(i));
    m2 = h * f(t(i) + (h/2), w(i) + (m1/2));
    m3 = h * f(t(i) + (h/2), w(i) + (m2/2));
    m4 = h * f(t(i+1), w(i) + m3); 
    w(i+1) = w(i) + (m1 + 2*m2 + 2*m3 + m4)/6;
  endfor

  %Algoritma utama
  for i = 4:n
    t(i+1) = t(i) + h;
  
  % Adams Bashfourth
    k1 = f(t(i), w(i));
    k2 = f(t(i-1), w(i-1));
    k3 = f(t(i-2), w(i-2));
    k4 = f(t(i-3), w(i-3));
    w(i+1) = w(i) + (h/24)*(55*k1-59*k2+37*k3-9*k4);
  
  % Adams Moulton
    k0 = f(t(i+1), w(i+1));
    w(i+1) = w(i) + (h/24)*(9*k0 + 19*k1 - 5*k2 + k3);
  endfor 
  