function [ x,iters,err ] = gauss_seidel( A,b,x,tol )
%Eric Matthews
%Gauss-Seidel Iterative Method
%HW 4 - NE 155 - March 2, 2015

iters = 0;
n = length(b);

L = tril(A);
U = A - L;

lastx = x;
x = inv(L)*(b - (U*x));
while (norm(x-lastx)/norm(x)) > tol
    lastx = x;
    x = inv(L)*(b - (U*x));
    iters = iters + 1;
end
err = (norm(x-lastx)/norm(x));

end