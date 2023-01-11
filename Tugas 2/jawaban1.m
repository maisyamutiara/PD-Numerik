clc;
clear all;
close all;
format long; %membuat 10 angka di belakang koma

function [xi, w1i, w2i] = linshoot(p, q, r, a, b, n, alpha, beta)
  %Inisialisasi
  h = (b-a)/n; %stepsize
  u = [alpha ; 0]; %nilai u(1,0) dan u(2,0)
  v = [0 ; 1];  %nilai v(1,0) dan v(2,0)
  xi = w1i = w2i = [];

  %Metode Runge-Kutta
  for i = 1:n
    x = a + (i-1)*h;

    k11 = h*u(2, i);
    k12 = h*(p(x)*u(2,i) + q(x)*u(1,i) + r(x));

    k21 = h*(u(2,i) + (k12/2));
    k22 = h*(p(x+(h/2))*(u(2,i)+(k12/2)) + q(x+(h/2))*(u(1,i)+(k11/2)) + r(x+(h/2)));

    k31 = h*(u(2,i) + (k22/2));
    k32 = h*(p(x+(h/2))*(u(2,i)+(k22/2)) + q(x+(h/2))*(u(1,i)+(k21/2)) + r(x+(h/2)));

    k41 = h*(u(2,i) + k32);
    k42 = h*(p(x+h)*(u(2,i)+k32) + q(x+h)*(u(1,i)+k31) + r(x+h));

    u(1, i+1) = u(1,i) + ((k11 + 2*k21 + 2*k31 + k41)/6);
    u(2, i+1) = u(2,i) + ((k12 + 2*k22 + 2*k32 + k42)/6);

    kp11 = h*v(2,i);
    kp12 = h*(p(x)*v(2,i) + q(x)*v(1,i));

    kp21 = h*(v(2,i) + (kp12/2));
    kp22 = h*(p(x+(h/2))*(v(2,i)+(kp12/2)) + q(x+(h/2))*(v(1,i)+(kp11/2)));

    kp31 = h*(v(2,i) + (kp22/2));
    kp32 = h*(p(x+(h/2))*(v(2,i)+(kp22/2)) + q(x+(h/2))*(v(1,i)+(kp21/2)));

    kp41 = h*(v(2,i) + kp32);
    kp42 = h*(p(x+h)*(v(2,i)+kp32) + q(x+h)*(v(1,i)+kp31));

    v(1, i+1) = v(1,i) + (kp11 + 2*kp21 + 2*kp31 + kp41)/6;
    v(2, i+1) = v(2,i) + (kp12 + 2*kp22 + 2*kp32 + kp42)/6;
  endfor

  w = [alpha; ((beta-u(1,(n+1)))/v(1,(n+1)))];
  xi(1) = a;
  w1i(1) = w(1,1);
  w2i(1) = w(2,1);

  %Solusi aproksimasi
  for i = 2:n+1
    W1 = u(1,i) + w(2,1)*v(1,i);
    W2 = u(2,i) + w(2,1)*v(2,i);
    x = a + (i-1)*h;
    xi(i) = x;
    w1i(i) = W1;
    w2i(i) = W2;
  endfor
endfunction

p = @(x) (x-1); %fungsi y'
q = @(x) (3*x^(-2)); %fungsi y
r = @(x) (x^(-1)*log(x)-1); %fungsi konstanta
a = 1; %batas kiri
b = 2; %batas kanan
alpha = 0; %nilai awal y(a)
beta = 0; %nilai batas y(b)
n = 10; %banyaknya meshpoint

[xi, w1i, w2i] = linshoot(p, q, r, a, b, n, alpha, beta);

[xi', w1i'] %print tabel nilai x dan solusinya

%visualisasi solusi aproksimasi
 hold on; %menampilkan beberapa plot dalam 1 gambar
 scatter(xi, w1i, 'b'); %scatter plot aproksimasi solusi
 legend('Aproksimasi'); %menambahkan legenda plot
 legend("location", "northwest");
 title('Linear Shooting Method'); %judul plot
