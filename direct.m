%Direct Solver - HW6 #6
%Eric Matthews
%April 6, 2015

%Matrix definition
a = 4;
h = 0.1;
Ea = 0.7;
D = 1;
%S0 = 8;
vEf = 0.6;

n = (2*a)/h + 1;

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

F = zeros(n);
for i = 1:n
    F(i,i) = vEf;
end

B = inv(A)*F;

[v,d] = eig(B);

k = max(max(d))

[row,col] = ind2sub(size(d),find(k==d));

phi1 = v(1:end,col);

%Power Iteration Algorithm
%1 Guess initials and Normalize
phi = ones(n,1);
k = 1;

phi = phi/norm(phi);

%2 (A is computed above)

%3 
Q = zeros(n,1);
for i = 1:n
    Q(i) = vEf*phi(i);
end

%4
E1 = 10^-4;
E2 = 10^-3;
conv_k = Inf;
conv_phi = Inf;
m=0;
while (conv_k > E1) && any(conv_phi > E2)
    m = m+1;
    %a
    phi_last = phi;
    D = diag(diag(A));
    L = tril(A) - D;
    U = A - D - L;
    w = 1.2;
    phi = inv(D + L.*w) * ((Q.*(1/k)).*w - (U.*w + D.*(w-1))*phi);
    
    %b
    Q_last = Q;
    Q = phi.*vEf;
    
    %c
    k_last = k;
    k = k*(sum(Q)/sum(Q_last));
    
    %d
    conv_k = abs((k-k_last)/k);
    conv_phi = zeros(n,1);
    for i = 1:n
        conv_phi(i) = abs((phi(i)-phi_last(i))/phi(i));
    end
end
phi2 = phi/norm(phi);

plot(-a:h:a,phi1)
set(gca,'FontSize',16)
xlabel('Position (cm)','FontSize',18)
ylabel('Flux (n/(cm^{3} s))','FontSize',18)

rel_err = zeros(n,1);
for i = 1:n
    rel_err(i) = abs((phi1(i) - phi2(i))/phi1(i));
end
norm(rel_err)

abs_err = abs(phi1 - phi2);
norm(abs_err)