clc;
clear all;
close all;

function [x,w] = linfd(p, q, r, aa, bb, alpha, beta, n)
    h = (bb - aa)/(n + 1); %stepsize
    x = aa + h; %interior mesh point

    %membuat matriks tridiagonal n*n dari sistem persamaan Ax = d
    a = b = c = d = []; %inisialisasi matriks kosong
    a(1) = 2 + (h^2) * q(x); %baris ke-1 kolom ke-1 matriks A
    b(1) = -1 + (h/2) * p(x); %baris ke-1 kolom ke-2 matriks A
    d(1) = -(h^2) * r(x) + (1+(h/2)*p(x)) * alpha; %baris ke-1 matriks d

    for i = 2:n-1
        x = aa + i*h;
        a(i) = 2 + (h^2) * q(x);
        b(i) = -1 + (h/2) * p(x);
        c(i) = -1 - (h/2) * p(x);
        d(i) = -(h^2) * r(x);
    endfor

    x = bb - h;
    a(n) = 2 + (h^2) * q(x); %baris ke-n kolom ke-n matriks A
    c(n) = -1 - (h/2) * p(x); %baris ke-n kolom ke-n-1 matriks A
    d(n) = -(h^2) * r(x) + (1 - (h/2)*p(x))*beta; %baris ke-n matriks d

    %Faktorisasi Crout untuk menyelesaikan tridiagonal sistem linear
    l = u = z = [];
    l(1) = a(1);
    u(1) = b(1) / a(1);
    z(1) = d(1) / l(1);

    %menyelesaikan masalah sistem linear pada matriks Lz = b
    for i = 2:n-1
        l(i) = a(i)-c(i)*u(i-1); %baris ke-i dari L
        u(i) = b(i)/l(i); %kolom ke i+1 dari U
        z(i) = (d(i)-c(i)*z(i-1))/l(i);
    endfor

    l(n) = a(n) - c(n)*u(n-1); %baris ke-n dari L
    z(n) = (d(n) - c(n)*z(n-1))/l(n);

    %menyelesaikan masalah sistem linear pada matriks Uw = z
    w = [];
    w(1) = alpha;
    w(n+2) = beta;
    w(n+1) = z(n);

    for i = 2:n
        w(n+2-i) = z(n+1-i)-u(n+1-i)*w(n+3-i);
    endfor

    %mencari nilai interior meshpoint dan menginputnya ke matriks
    x = [];
    for i = 1:n+2
        x(i) = aa + (i-1)*h;
    endfor
endfunction

p = @(x) (2); %koefisien dari y' yaitu p(x)y'
q = @(x) (-1); %koefisien dari y yaitu q(x)y
r = @(x) (x*exp(x)-x); %konstanta
aa = 0; %batas kiri
bb = 2; %batas kanan
alpha = 0; %nilai awal y(aa)
beta = -4; %nilai batas y(bb)
n = 9; %banyaknya mesh point

[x,w] = linfd(p, q, r, aa, bb, alpha, beta, n);

%pendefinisian fungsi eksak
sol = @(x) ((1/6).*x.^3.*exp(x) - (5/3).*x.*exp(x) + 2.*exp(x) - x - 2);

format long %membuat 10 angka di belakang koma
[x', w', sol(x)', abs(sol(x)-w)'] %print tabel nilai x, aproksimasi, solusi eksak, dan error

%visualisasi solusi aproksimasi
hold on; %menampilkan beberapa plot dalam 1 gambar
fplot(sol, [0,2], 'k'); %plot fungsi eksak
scatter(x, w, 'r'); %scatter plot aproksimasi solusi
legend('Eksak', 'Aproksimasi', 'location', 'northwest'); %menambahkan legenda plot
title('Finite Difference Linear Methods'); %judul plot
