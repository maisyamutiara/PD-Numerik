clc;
clear all;
close all;
format long;

function [x, t, u] = wave(c, phi, psi, lb, rb, xb, xu, tb, tu, dx, dt)
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  t = tb:dt:tu; %partisi tb sampai tu sebesar stepsize dt
  u = [];
  nt = length(t); %menyimpan panjang Nt
  nx = length(x); %menyimpan panjang Nx
  S = (c^2 * dt^2) / dx^2;

  %nilai awal u_j^0
  for j = 1:nx
    u(j, 1) = phi(x(j));
  endfor

  %nilai awal u_j^1
  for j = 2:nx-1
    u(j, 2) = (S/2) * (phi(x(j+1)) + phi(x(j-1))) + (1-S) * phi(x(j)) + dt * psi(x(j));
  endfor

  %syarat batas
  for n = 2:nt
    u(1, n) = lb(t(n));
    u(nx, n) = rb(t(n));
  endfor

  %persamaan beda u_j^(n+1)
  for n = 2:nt-1
    for j = 2:nx-1
      u(j, n+1) = S * (u(j+1, n) + u(j-1, n)) + 2 * (1-S) * u(j, n) - u(j, n-1);
    endfor
  endfor
endfunction

c = 1; %nilai c pada utt - c^2 uxx = 0
phi = @(x) sin(2*pi*x); %nilai awal u(x,0)
psi = @(x) 2*pi*sin(2*pi*x); %nilai awal ut(x,0)
lb = rb = @(t) 0; %syarat batas u(0,t) = u(1,t)
xb = 0; %batas bawah variabel x
xu = 1; %batas atas variabel x
tb = 0; %batas bawah variabel t
tu = 1; %batas atas variabel t
dx = 0.1; %stepsize h
dt = 0.1; %stepsize k

%solusi aproksimasi
[x, t, u] = wave(c, phi, psi, lb, rb, xb, xu, tb, tu, dx, dt);
[x'] %print nilai x
[t'] %print nilai t
[u'] %print solusi aproksimasi

%solusi eksak
sol = @(x, t) (sin(2*pi*x)*(cos(2*pi*t)+sin(2*pi*t)));

for j = 1:length(x)
  for n = 1:length(t)
    y(j, n) = sol(x(j), t(n));
  endfor
endfor

w1 = @(x) (sin(2.*pi.*x).*(cos(2.*pi)+sin(2.*pi))); %solusi eksak saat t = 1
[w1(x)'] %print solusi eksak

figure(1);
mesh(x, t, u');
xlabel("x");
ylabel("t");
zlabel("u");
title("Solusi aproksimasi");

figure(2);
mesh(x, t, y');
xlabel("x");
ylabel("t");
zlabel("u");
title("Solusi eksak");

figure(3);
hold on;
fplot(w1, [0, 1], 'k');
plot(x, u(:, length(t)), 'r');
legend('Eksak', 'Aproksimasi');
title("Perbandingan solusi eksak dan solusi aproksimasi pada t=1");
