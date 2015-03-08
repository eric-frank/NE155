function [ x,iters,err ] = SOR( A,b,w,x,tol )
%Eric Matthews
%SOR Iterative Method
%HW 4 - NE 155 - March 2, 2015

iters = 0;
n = length(b);

D = diag(diag(A));
L = tril(A) - D;
U = A - D - L;

lastx = x;
x = inv(D + L.*w) * (b.*w - (U.*w + D.*(w-1))*x);
while (norm(x-lastx)/norm(x)) > tol
    lastx = x;
    x = inv(D + L.*w) * (b.*w - (U.*w + D.*(w-1))*x);
    iters = iters + 1;
end
err = (norm(x-lastx)/norm(x));

end

