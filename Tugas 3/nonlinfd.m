clc;
clear all;
close all;

function [xi, w] = nonlinfd(f, fy, fyp, aa, bb, n, alpha, beta, m, tol)
  h = (bb - aa)/(n+1); %stepsize
  xi = w = [];
  w(1) = alpha;
  w(n+2) = beta;

  for i = 2:n+1
    w(i) = alpha + (i-1)*((beta-alpha)/(bb-aa))*h;
  endfor

  k = 1;
  while k <= m
    x = aa + h;
    t = (w(3)-alpha) / (2*h);
    a = b = c = d = l = u = z = [];
    a(1) = 2 + (h^2) * fy(x, w(2), t);
    b(1) = -1 + (h/2) * fyp(x, w(2), t);
    d(1) = -(2*w(2) - w(3) - alpha + (h^2)*f(x, w(2), t));

    for i = 2:n-1
      x = aa + i*h;
      t = (w(i+2) - w(i)) / (2*h);
      a(i) = 2 + (h^2) * fy(x, w(i+1), t);
      b(i) = -1 + (h/2) * fyp(x, w(i+1), t);
      c(i) = -1 - (h/2) * fyp(x, w(i+1), t);
      d(i) = -(2*w(i+1) - w(i+2) - w(i) + (h^2)*f(x, w(i+1), t));
    endfor

    x = bb - h;
    t = (beta-w(n))/(2*h);
    a(n) = 2 + (h^2) * fy(x, w(n+1), t);
    c(n) = -1 - (h/2) * fyp(x, w(n+1), t);
    d(n) = -(2*w(n+1) - w(n) - beta + (h^2)*f(x, w(n+1), t));

    %Faktorisasi Crout untuk menyelesaikan tridiagonal sistem linear
    l(1) = a(1);
    u(1) = b(1)/a(1);
    z(1) = d(1)/l(1);

    for i = 2:n-1
      l(i) = a(i) - c(i) * u(i-1);
      u(i) = b(i)/l(i);
      z(i) = (d(i) - c(i)*z(i-1))/l(i);
    endfor

    l(n) = a(n) - c(n) * u(n-1);
    z(n) = (d(n) - c(n)*z(n-1))/l(n);
    v(n) = z(n);
    w(n+1) = w(n+1)+v(n);

    for i = n-1:-1:1
      v(i) = z(i) - u(i)*v(i+1);
      w(i+1) = w(i+1) + v(i);
    endfor

    if norm(v) <= tol
      for i = 1:n+2
        xi(i) = aa + (i-1)*h;
      endfor

      return
    endif

    k = k+1;
  endwhile

  disp('max iteration')
endfunction

f = @(x, y, yp) (32 + 2*x^3 - y*yp)/8; %fungsi y"
fy = @(x, y, yp) -18*y/x^5; %turunan parsial fungsi y" terhadap y
fyp = @(x, y, yp) 2*yp/x^3; %turunan parsial fungsi y" terhadap y'
aa = 1; %batas kiri
bb = 2; %batas kanan
alpha = 0; %nilai awal y(aa)
beta = log(256); %nilai batas y(bb)
n = 19; %banyaknya meshpoint
m = 20; %banyaknya maksimum iterasi
tol = 10^-4; %batas toleransi error

[xi, w] = nonlinfd(f, fy, fyp, aa, bb, n, alpha, beta, m, tol);

sol = @(x) (x.^3.*log(x)); %pendefinisian fungsi eksak

format long %membuat 10 angka di belakang koma
[xi', w', sol(xi)', abs(sol(xi)-w)'] %print tabel nilai x, aproksimasi, solusi eksak, dan error

%visualisasi solusi aproksimasi
hold on; %menampilkan beberapa plot dalam 1 gambar
fplot(sol, [1, 2], 'k'); %plot fungsi eksak
scatter(xi, w, 'r'); %scatter plot aproksimasi solusi
legend('Eksak', 'Aproksimasi', 'location', 'northwest'); %menambahkan legenda plot
title('Finite Difference Nonlinear Methods'); %judul plot
