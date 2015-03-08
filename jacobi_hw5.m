function [ x,iters,err ] = jacobi( A,b,guess,tol )
%Eric Matthews
%Jacobi Iterative Method
%HW 4 - NE 155 - March 2, 2015

iters = 1;
n = length(b);

D = diag(diag(A));
R = A-D;

Pj = eye(n) - (inv(D)*A);

if max(abs(eig(Pj))) < 1
    lastx = guess;
    x = ((inv(D))*(D-A)*guess) + (inv(D) * b);
    while (norm(x-lastx)/norm(x)) > tol
        lastx = x;
        x = ((inv(D))*(D-A)*x) + (inv(D) * b);
        iters = iters + 1;
    end
    err = (norm(x-lastx)/norm(x));
else
    error('NO CONVERGENCE')
end

end
