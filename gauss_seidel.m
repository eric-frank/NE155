function [ x,iters ] = gauss_seidel( A,b,x,tol )
%Eric Matthews
%Gauss-Seidel Iterative Method
%HW 4 - NE 155 - March 2, 2015

iters = 0;
n = length(b);

L = tril(A);
U = A - L;

while norm(A*x - b) > tol
    x = inv(L)*(b - (U*x));
    iters = iters + 1;
end

end