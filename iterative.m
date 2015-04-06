%Iterative Methods HW6 - #4
%Eric Matthews
%April 6, 2015

%Matrix definition
a = 4;
h = 0.1;
Ea = 0.2;
D = 1;
S0 = 8;

n = (2*a)/h + 1;

S = [];
for i = 1:n
    S(i) = S0;
end


A = zeros(n);

A(1,1) = ((2*D)/(h^2)) + Ea;
A(1,2) = (-1*D)/(h^2);

for i = 2:n-1
    A(i,i) = ((2*D)/(h^2)) + Ea;
    A(i,i-1) = (-1*D)/(h^2);
    A(i,i+1) = (-1*D)/(h^2);
end

A(n,n-1) = (-1*D)/(h^2);
A(n,n) = ((2*D)/(h^2)) + Ea;



%Jacobi
guess = zeros(n,1);
[phi,iters] = jacobi(A,S',guess,10^-3);
iters

plot(-a:h:a,phi)
set(gca,'FontSize',16)
xlabel('Position (cm)','FontSize',18)
ylabel('Flux (n/(cm^{3} s))','FontSize',18)

%pause

[phi,iters] = jacobi(A,S',guess,10^-5);
iters

plot(-a:h:a,phi)
set(gca,'FontSize',16)
xlabel('Position (cm)','FontSize',18)
ylabel('Flux (n/(cm^{3} s))','FontSize',18)



%Gauss-Seidel
[phi,iters] = gauss_seidel(A,S',guess,10^-3);
iters

plot(-a:h:a,phi)
set(gca,'FontSize',16)
xlabel('Position (cm)','FontSize',18)
ylabel('Flux (n/(cm^{3} s))','FontSize',18)

%pause

[phi,iters] = gauss_seidel(A,S',guess,10^-5);
iters

plot(-a:h:a,phi)
set(gca,'FontSize',16)
xlabel('Position (cm)','FontSize',18)
ylabel('Flux (n/(cm^{3} s))','FontSize',18)



%SOR
[phi,iters] = SOR(A,S',1.2,guess,10^-3);
iters

plot(-a:h:a,phi)
set(gca,'FontSize',16)
xlabel('Position (cm)','FontSize',18)
ylabel('Flux (n/(cm^{3} s))','FontSize',18)

%pause

[phi,iters] = SOR(A,S',1.2,guess,10^-5);
iters

plot(-a:h:a,phi)
set(gca,'FontSize',16)
xlabel('Position (cm)','FontSize',18)
ylabel('Flux (n/(cm^{3} s))','FontSize',18)