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
 hold on;
 scatter(xi, w1i, 'b');
 legend('Aproksimasi');
 legend("location", "northwest");
 title('Linear Shooting Method');
