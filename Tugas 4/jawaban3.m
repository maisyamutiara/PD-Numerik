clc;
clear all;
close all;
format long;

function [x, t, u] = explicitheat(d, f, lb, rb, xb, xu, tb, tu, dx, dt)
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  t = tb:dt:tu; %partisi tb sampai tu sebesar stepsize dt
  u = [];
  nt = length(t); %menyimpan panjang Nt
  nx = length(x); %menyimpan panjang Nx
  S = (d * dt) / dx^2;

  %nilai awal u_j^0
  for j = 1:nx
    u(j, 1) = f(x(j));
  endfor

  %syarat batas
  for n = 2:nt
    u(1, n) = lb(t(n));
    u(nx, n) = rb(t(n));
  endfor

  %persamaan beda hingga u_j^(n+1)
  for n = 1:nt-1
    for j = 2:nx-1
      u(j, n+1) = (1 - 2*S) * u(j, n) + S * (u(j+1, n) + u(j-1, n));
    endfor
  endfor
endfunction

d = 1; %nilai d pada ut - d uxx = 0
f = @(x) sin(2*pi*x); %nilai awal u(x,0)
lb = rb = @(t) 0; %syarat batas u(0,t) = u(1,t)
xb = 0; %batas bawah variabel x
xu = 2; %batas atas variabel x
tb = 0; %batas bawah variabel t
tu = 1; %batas atas variabel t
dx = 0.2; %stepsize h
dt = 0.1; %stepsize k

[x, t, u] = explicitheat(d, f, lb, rb, xb, xu, tb, tu, dx, dt);
[x'] %print nilai x
[t'] %print nilai t
[u'] %print solusi aproksimasi

%solusi eksak
sol = @(x, t) (exp(-4*pi^2*t)*sin(2*pi.*x))

for j = 1:length(x)
  for n = 1:length(t)
    y(j, n) = sol(x(j), t(n));
    [y(j, n)']
  endfor
endfor

[sol(x, 1)'] %print solusi eksak

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
