function [ x,iters ] = SOR( A,b,w,x,tol )
%Eric Matthews
%SOR Iterative Method
%HW 4 - NE 155 - March 2, 2015

iters = 0;
n = length(b);

D = diag(diag(A));
L = tril(A) - D;
U = A - D - L;

while norm(A*x -b) > tol
    x = inv(D + L.*w) * (b.*w - (U.*w + D.*(w-1))*x);
    iters = iters + 1;
end

end

