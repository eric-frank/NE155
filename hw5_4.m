n = 5;

A = zeros(n);
b = [];

for i = 1:n
    A(i,i) = 4;
    b(i) = 100;
end
b = b';

for i = 1:(n-1)
    A(i+1,i) = -1;
    A(i,i+1) = -1;
end

guess = [];
for i = 1:n
    guess(i) = 0;
end
guess = guess';

inv(A)*b

tol = 10^-6;

[x,iters,err] = jacobi_hw5(A,b,guess,tol);
err
[x,iters,err] = gauss_seidel_hw5(A,b,guess,tol);
err
[x,iters,err] = SOR_hw5(A,b,1.1,guess,tol);
err

[x,iters,err] = jacobi(A,b,guess,tol);
err
[x,iters,err] = gauss_seidel(A,b,guess,tol);
err
[x,iters,err] = SOR(A,b,1.1,guess,tol);
err