clc;
clear all;
close all;
format long g;

function [x, t, u] = courant(d, f, lb, xb, xu, tb, tu, dx, dt)
  t = tb:dt:tu; %partisi tb sampai tu sebesar stepsize dt
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  u = [];

  %nilai awal u_j^0
  for j = 1:length(x)
    u(j, 1) = f(x(j));
  end

  %syarat batas
  for n = 1:length(t)
    u(1, n) = lb(t(n));
  end

  %persamaan beda hingaa u_j^(n+1)
  c = d * dt / dx;
  for n = 1:length(t)-1
    for j = 2:length(x)
      u(j, n+1) = (1-c) * u(j, n) + c * u(j-1, n);
    end
  end
end

function [x, t, u] = richardson(d, f, lb, rb, xb, xu, tb, tu, dx, dt)
  t = tb:dt:tu; %partisi tb sampai tu sebesar stepsize dt
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  nt = length(t); %menyimpan panjang Nt
  nx = length(x); %menyimpan panjang Nx
  u = [];

  %nilai awal u_j^0
  for j = 1:nx
    u(j, 1) = f(x(j));
  end

  %syarat batas
  for n = 1:nt
    u(1, n) = lb(t(n));
    u(nx, n) = rb(t(n));
  end

  %persamaan beda hingga u_j^(n+1)
  c = (d*dt) / (2*dx);
  for n = 1:nt-1
    for j = 2:nx-1
      u(j, n+1) = u(j, n) - (c * (u(j+1, n) - u(j-1, n)));
    end
  end
end

function [x, t, u] = lax(d, f, lb, rb, xb, xu, tb, tu, dx, dt)
  t = tb:dt:tu; %partisi tb sampai tu sebesar stepsize dt
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  nt = length(t); %menyimpan panjang Nt
  nx = length(x); %menyimpan panjang Nx
  u = [];

  %nilai awal u_j^0
  for j = 1:nx
    u(j, 1) = f(x(j));
  end

  %syarat batas
  for n = 1:nt
    u(1, n) = lb(t(n));
    u(nx, n) = rb(t(n));
  end

  %persamaan beda hingga u_j^(n+1)
  c = (d*dt) / (2*dx);
  for n = 1:nt-1
    for j = 2:nx-1
      u(j, n+1) = ((u(j+1, n) + u(j-1, n)) / 2) - (c*(u(j+1, n) - u(j-1, n)));
    end
  end
end

d = 2; %nilai d pada ut - d uxx = 0
f = @(x) exp(-x^2); %nilai awal u(x,0)
lb = @(t) exp(-4*t^2); %syarat batas u(0,t)
rb = @(t) exp(-1*((2-2*t)^2)); %syarat batas u(2,t)
xb = 0; %batas bawah variabel x
xu = 2; %batas atas variabel x
tb = 0; %batas bawah variabel t
tu = 1; %batas atas variabel t
dx = 0.2; %stepsize h
dt = 0.1; %stepsize k

%solusi aproksimasi
[x, t, uc] = courant(d, f, lb, xb, xu, tb, tu, dx, dt);
[x, t, ur] = richardson(d, f, lb, rb, xb, xu, tb, tu, dx, dt);
[x, t, ul] = lax(d, f, lb, rb, xb, xu, tb, tu, dx, dt);
[x'] %print nilai x
[t'] %print nilai t
[uc'] %print solusi aproksimasi metode courant
[ur'] %print solusi aproksimasi metode richardson
[ul'] %print solusi aproksimasi metode lax

%solusi eksak
sol = @(x, t) exp(-1*((x-2*t)^2));

for j = 1:length(x)
  for n = 1:length(t)
    y(j, n) = sol(x(j), t(n));
  endfor
endfor

w1 = @(x) exp(-1.*((x-2).^2)); %solusi eksak saat t = 1
[w1(x)'] %print solusi eksak

figure(1);
hold on;
mesh(x, t, uc');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Aproksimasi Metode Courant-Isaacson-Rees");

figure(2);
hold on;
mesh(x, t, ur');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Aproksimasi Metode Richardson");

figure(3);
hold on;
mesh(x, t, ul');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Aproksimasi Metode Lax");

figure(4);
hold on;
mesh(x, t, y');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Eksak")

figure(5);
hold on;
fplot(w1, [0, 2], 'k');
plot(x, uc(:, length(t)), 'r');
xlabel('x');
ylabel('u');
legend("Eksak", "Metode Courant");
title("Eksak vs Courant di penampang pada t=1");

figure(6);
hold on;
fplot(w1, [0, 2], 'k');
plot(x, ur(:, length(t)), 'r');
xlabel('x');
ylabel('u');
legend("Eksak", "Metode Richardson");
title("Eksak vs Richardson di penampang pada t=1");

figure(7);
hold on;
fplot(w1, [0, 2], 'k');
plot(x, ul(:, length(t)), 'r');
xlabel('x');
ylabel('u');
legend("Eksak", "Metode Lax");
title("Eksak vs Lax di penampang pada t=1");
