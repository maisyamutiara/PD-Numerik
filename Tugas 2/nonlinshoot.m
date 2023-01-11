function [xi, w1i, w2i] = nonlinshoot(f, fy, fyp, a, b, n, alpha, beta, m, tol)
  h = (b-a)/n;
  k = 1;
  tk = (beta-alpha)/(b-a);
  xi = w1i = w2i = [];
  while k <= m
    w = [alpha; tk];
    u = [0,1];
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

    if abs(w(1, n+1)-beta) <= tol
      for i = 1:n+1
        x = a +(i-1)*h;
        xi(i) = x;
        w1i(i) = w(1,i);
        w2i(i) = w(2,i);
      endfor
      return
    endif
    tk = tk - ((w(1,n+1)-beta)/u(1));
    k = k+1;
  endwhile
  disp('max iteration');
endfunction
