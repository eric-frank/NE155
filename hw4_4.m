%part a

n = 100;

A = zeros(n);
b = [];

A(1,1) = 2;
A(1,2) = -1;
b(1) = 0;

for i = 2:(n-1)
    A(i,i-1) = -1;
    A(i,i) = 2;
    A(i,i+1) = -1;
    b(i) = i-1;
end

A(n,n-1) = -1;
A(n,n) = 2;
b(n) = 99;

%part b

condition = norm(inv(A)) * norm(A)

%part c

x1 = inv(A) * b'

%part d

x2 = A\b'

%plotting

plot(b,x1)
hold on
plot(b,x2,'go')
xlabel('b','FontSize',18)
ylabel('x','FontSize',18)
legend('Explicit','Backslash Operator')



