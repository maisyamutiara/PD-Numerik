clc;
clear all;
close all;
format long;

function [x, y, u] = LaplaceJacobi(lb, rb, ub, db, xb, xu, yb, yu, dx, dy, N)
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  y = yb:dy:yu; %partisi yb sampai tu sebesar stepsize dy
  nx = length(x); %menyimpan panjang Nx
  ny = length(y); %menyimpan panjang Ny

  r = 0.5*dy^2/(dy^2+dx^2);
  s = 0.5*dx^2/(dy^2+dx^2);

  %nilai awal u_j^0
  for i = 1:nx
    u(i, 1) = db(x(i));
    u(i, ny) = ub(x(i));
  endfor

  %syarat batas
  for j = 2:ny
    u(1, j) = lb(y(j));
    u(nx, j) = rb(y(j));
  endfor

  %Algoritma Jacobi
  u2 = u;

  for n = 1:N
    for i = 2:nx-1
      for j = 2:ny-1
        u2(i, j) = r * (u(i+1, j) + u(i-1, j)) + s * (u(i, j+1) + u(i, j-1));
      endfor
    endfor
    for i = 2:nx-1
      for j = 2:ny-1
        u(i, j) = u2(i, j);
      endfor
    endfor
  endfor
endfunction

function [x, y, u] = LaplaceGS(lb, rb, ub, db, xb, xu, yb, yu, dx, dy, N)
  x = xb:dx:xu; %partisi xb sampai xu sebesar stepsize dx
  y = yb:dy:yu; %partisi yb sampai tu sebesar stepsize dy
  nx = length(x); %menyimpan panjang Nx
  ny = length(y); %menyimpan panjang Ny

  r = 0.5*dy^2/(dy^2+dx^2);
  s = 0.5*dx^2/(dy^2+dx^2);

  %nilai awal u_j^0
  for i = 1:nx
    u(i, 1) = db(x(i));
    u(i, ny) = ub(x(i));
  endfor

  %syarat batas
  for j = 2:ny
    u(1, j) = lb(y(j));
    u(nx, j) = rb(y(j));
  endfor

  %persamaan beda u_j^(n+1)
  for n = 1:N
    for i = 2:nx-1
      for j = 2:ny-1
        u(i, j) = r * (u(i+1, j) + u(i-1, j)) + s * (u(i, j+1) + u(i, j-1));
      endfor
    endfor
  endfor
endfunction

db = @(x) 2*log(x); %nilai awal u(x,0)
ub = @(x) log((x^2)+1); %nilai awal u(x,1)
lb = @(y) log((y^2)+1); %syarat batas u(1,y)
rb = @(y) log((y^2)+4); %syarat batas u(2,y)
xb = 1; %batas bawah variabel x
xu = 2; %batas atas variabel x
yb = 0; %batas bawah variabel y
yu = 1; %batas atas variabel y
dx = 1/3; %stepsize h
dy = 1/3; %stepsize k
N = 20; %Iterasi maksimum

[x, y, ujc] = LaplaceJacobi(lb, rb, ub, db, xb, xu, yb, yu, dx, dy, N);
[x, y, ugs] = LaplaceGS(lb, rb, ub, db, xb, xu, yb, yu, dx, dy, N);
[x'] %print nilai x
[y'] %print nilai y
[ujc'] %print solusi aproksimasi jacobi
[ugs'] %print solusi aproksimasi gauss-seidel

%Solusi eksak
sol = @(x, y) (log((x^2)+(y^2)))

for m = 1:length(x)
  for n = 1:length(y)
    w(m, n) = sol(x(m), y(n));
    [w(m, n)'] %print solusi eksak
  endfor
endfor

figure(1);
mesh(x, y, ujc');
xlabel("x");
ylabel("y");
zlabel("u");
title("Solusi Aproksimasi Metode Iterasi Jacobi");

figure(2);
mesh(x, y, ugs');
xlabel("x");
ylabel("y");
zlabel("u");
title("Solusi Aproksimasi Metode Iterasi Gauss-Seidel");

figure(3);
mesh(x, y, w');
xlabel("x");
ylabel("y");
zlabel("u");
title("Solusi eksak");
