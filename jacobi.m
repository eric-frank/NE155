function [ x,iters ] = jacobi( A,b,guess,tol )
%Eric Matthews
%Jacobi Iterative Method
%HW 4 - NE 155 - March 2, 2015

iters = 1;
n = length(b);

D = diag(diag(A));
R = A-D;

Pj = eye(n) - (inv(D)*A);

if max(abs(eig(Pj))) < 1
    x = ((inv(D))*(D-A)*guess) + (inv(D) * b);
    while norm(A*x -b) > tol
        x = ((inv(D))*(D-A)*x) + (inv(D) * b);
        iters = iters + 1;
    end
else
    error('NO CONVERGENCE')
end

end
