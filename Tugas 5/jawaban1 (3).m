clc;
clear all;
close all;
format long;

function [x, t, w] = ImplicitHeat(d, f, lb, rb, xb, xu, tb, tu, dx, dt)
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  t = tb:dt:tu; %partisi tb sampai tu sebesar stepsize dt
  nx = length(x); %menyimpan panjang Nx
  nt = length(t); %menyimpan panjang Nt

  # Nilai lambda
  lambd = (d * dt) / (dx^2);

  # Nilai awal dan syarat batas
  for i = 1:nx
    w(i, 1) = f(x(i));
  endfor

  for j = 2:nt
    w(1, j) = lb(t(j));
    w(nx, j) = rb(t(j));
  endfor

  # Penyelesaian SPL dengan faktorisasi Crout
  l(2) = 1 + 2*lambd;
  u(2) = -lambd / l(2);
  for i = 3:nx-2
    l(i) = 1 + 2*lambd + lambd*u(i-1);
    u(i) = -lambd / l(i);
  endfor
  l(nx-1) = 1 + 2*lambd + lambd*u(nx-2);
  for j = 2:nt
    z(2) = w(2, j-1) / l(2);
    for i = 3:nx-1
      z(i) = (w(i, j-1) + lambd*z(i-1)) / l(i);
    endfor
    w(nx-1, j) = z(nx-1);
    for i = nx-2:-1:2
      w(i, j) = z(i) - u(i)*w(i+1, j);
    endfor
  endfor
endfunction

function [x, t, w] = CrankNicolson(d, f, lb, rb, xb, xu, tb, tu, dx, dt)
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  t = tb:dt:tu; %partisi tb sampai tu sebesar stepsize dt
  nx = length(x); %menyimpan panjang Nx
  nt = length(t); %menyimpan panjang Nt

  # Nilai lambda
  lambd = (d * dt) / (dx^2);

  # Nilai awal dan syarat batas
  for i = 1:nx
    w(i, 1) = f(x(i));
  endfor

  for j = 2:nt
    w(1, j) = lb(t(j));
    w(nx, j) = rb(t(j));
  endfor

  # Penyelesaian SPL menggunakan faktorisasi Crout
  l(2) = 1 + lambd;
  u(2) = -lambd / (2*l(2));
  for i = 3:nx-2
    l(i) = 1 + lambd + (lambd*u(i-1))/2;
    u(i) = -lambd / (2*l(i));
  endfor
  l(nx-1) = 1 + lambd + (lambd*u(nx-2))/2;
  for j = 2:nt
    z(2) = ((1-lambd)*w(2, j-1) + (lambd/2)*w(3, j-1)) / l(2);
    for i = 3:nx-1
      z(i) = ((1-lambd)*w(i, j-1) + (lambd/2)*(w(i+1, j-1) + w(i-1, j-1) + z(i-1))) / l(i);
    endfor
    w(nx-1, j) = z(nx-1);
    for i = nx-2:-1:2
      w(i, j) = z(i) - u(i)*w(i+1, j);
    endfor
  endfor
endfunction

d = 4/pi^2; %nilai d pada ut - d uxx = 0
f = @(x) sin(pi*x/4)*(1+2*cos(pi*x/4)); %nilai awal u(x,0)
lb = rb = @(t) 0; %syarat batas u(0,t) = u(4,t)
xb = 0; %batas bawah variabel x
xu = 4; %batas atas variabel x
tb = 0; %batas bawah variabel t
tu = 1; %batas bawah variabel t
dx = 0.4; %stepsize h
dt = 0.1; %stepsize k

%solusi aproksimasi
[x, t, wbtcs] = ImplicitHeat(d, f, lb, rb, xb, xu, tb, tu, dx, dt);
[x, t, wcn] = CrankNicolson(d, f, lb, rb, xb, xu, tb, tu, dx, dt);
[x'] %print nilai x
[t'] %print nilai t
[wbtcs'] %print solusi aproksimasi BTCS
[wcn'] %print solusi aproksimasi Crank Nicolson

%solusi eksak
u = @(x, t) exp(-t)*sin(x*pi/2) + exp(-t/4)*sin(x*pi/4);
for i = 1:length(x)
  for j = 1:length(t)
    ufig(i, j) = u(x(i), t(j));
  endfor
endfor
u1 = @(x, t) exp(-1)*sin(x*pi./2) + exp(-1/4)*sin(x*pi/4); %solusi eksak saat t = 1
[u1(x, t)'] %print solusi eksak saat t = 1

figure(1);
mesh(x, t, ufig');
xlabel("x");
ylabel("t");
zlabel("u");
title("Solusi eksak");

figure(2);
mesh(x, t, wbtcs');
xlabel("x");
ylabel("t");
zlabel("u");
title("Solusi aproksimasi BTCS");

figure(3);
mesh(x, t, wcn');
xlabel("x");
ylabel("t");
zlabel("u");
title("Solusi aproksimasi Crank Nicolson");

figure(4);
hold on;
fplot(u1, [0, 4], 'k');
plot(x, wbtcs(:, length(t)), 'r');
plot(x, wcn(:, length(t)), 'c');
legend('Eksak', 'BTCS', 'Crank Nicolson');
title("Perbandingan solusi eksak dan solusi aproksimasi pada t=1");
