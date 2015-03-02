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

[x,iters] = jacobi(A,b,guess,10^-6)
[x,iters] = gauss_seidel(A,b,guess,10^-6)
[x,iters] = SOR(A,b,1.1,guess,10^-6)