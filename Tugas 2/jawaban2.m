clc;
clear all;
close all;
format long; %membuat 10 angka di belakang koma

function [xi, w1i, w2i] = nonlinshoot(f, fy, fyp, a, b, n, alpha, beta, m, tol)
  %inisialisasi
  h = (b-a)/n; %stepsize
  k = 1;
  tk = (beta-alpha)/(b-a); %guess
  xi = w1i = w2i = [];
  while k <= m
    w = [alpha; tk]; %nilai w(1,0) dan w(2,0)
    u = [0,1]; %nilai u(1) dan u(2)

    %Metode Runge-Kutta
    for i = 1:n
      x = a + (i-1)*h;
      k11 = h*w(2,i);
      k12 = h*f(x, w(1,i), w(2,i));

      k21 = h*(w(2,i) + (k12/2));
      k22 = h*f((x+(h/2)), (w(1,i)+(k11/2)), (w(2,i)+(k12/2)));

      k31 = h*(w(2,i) + (k22/2));
      k32 = h*f((x+(h/2)), (w(1,i)+(k21/2)), (w(2,i)+(k22/2)));

      k41 = h*(w(2,i) + k32);
      k42 = h*f((x+h), (w(1,i)+k31), (w(2,i)+k32));

      w(1, i+1) = w(1,i) + (k11 + 2*k21 + 2*k31 + k41)/6;
      w(2, i+1) = w(2,i) + (k12 + 2*k22 + 2*k32 + k42)/6;

      kp11 = h*u(2);
      kp12 = h*(fy(x, w(1,i), w(2,i))*u(1) + fyp(x, w(1,i), w(2,i))*u(2));

      kp21 = h*(u(2) + (kp12/2));
      kp22 = h*(fy((x+(h/2)), w(1,i), w(2,i))*u(1) + fyp((x+(h/2)), w(1,i), w(2,i))*(u(2) + (kp12/2)));

      kp31 = h*(u(2) + (kp22/2));
      kp32 = h*(fy((x+(h/2)), w(1,i), w(2,i))*(u(1) + (kp21/2)) + fyp((x+(h/2)), w(1,i), w(2,i))*(u(2) + (kp22/2)));

      kp41 = h*(u(2) + kp32);
      kp42 = h*(fy((x+h), w(1,i), w(2,i))*(u(1) + kp31) + fyp((x+h), w(1,i), w(2,i))*(u(2) + kp32));

      u(1) = u(1) + (kp11 + 2*kp21 + 2*kp31 + kp41)/6;
      u(2) = u(2) + (kp12 + 2*kp22 + 2*kp32 + kp42)/6;
    endfor

    %Periksa batas error
    if abs(w(1, n+1)-beta) <= tol
      for i = 1:n+1
        x = a +(i-1)*h;
        xi(i) = x;
        w1i(i) = w(1,i);
        w2i(i) = w(2,i);
      endfor
      return
    endif

    %Metode Newton untuk menghitung tk
    tk = tk - ((w(1,n+1)-beta)/u(1));
    k = k+1;
  endwhile
  disp('max iteration');
endfunction

f = @(x, y, yp) ((1/2)*(1-((yp)^2)- y*sin(x))); %fungsi y"
fy = @(x, y, yp) (-(sin(x))/2); %turunan parsial fungsi y" terhadap y
fyp = @(x, y, yp) (-yp); %turunan parsial fungsi y" terhadap y'
a = 0; %batas kiri
b = pi; %batas kanan
alpha = 2; %nilai awal y(a)
beta = 2; %nilai batas y(b)
n = 20; %banyaknya meshpoint
m = 10; %banyaknya maksimum iterasi
tol = 10^(-4); %batas toleransi error

[xi, w1i, w2i] = nonlinshoot(f, fy, fyp, a, b, n, alpha, beta, m, tol);

%pendefinisian fungsi eksak
sln = @(x) (2 + sin(x));
wex = [];
for i = 1:length(xi)
  wex(i) = sln(xi(i));
endfor

[xi', w1i', wex', abs(wex-w1i)']%print tabel nilai x, aproksimasi, solusi eksak, dan error

%visualisasi solusi aproksimasi
 hold on; %menampilkan beberapa plot dalam 1 gambar
 fplot(sln, [0,pi], 'k'); %plot fungsi eksak
 scatter(xi, w1i, 'b'); %scatter plot aproksimasi solusi
 legend('Eksak','Aproksimasi'); %menambahkan legenda plot
 legend("location", "northwest");
 title('Nonlinear Shooting Method'); %judul plot

