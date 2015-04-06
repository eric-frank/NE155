%Numerical to Analytical Comparison #HW6 #3
%Eric Matthews
%April 6, 2015

err = [];
spaces = [];
for h = [1,0.5,0.1,0.05,0.01]
    %Matrix definition
    a = 4;
    Ea = 0.2;
    D = 1;
    L = sqrt(D/Ea);
    S0 = 8;
    
    as = -4:h:4;
    n = length(as);
    spaces = [spaces,length(-a:h:a)];

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
    
    %Algorithm
    u = zeros(n,1);
    v = zeros(n,1);

    %1
    u(1) = A(1,1);
    v(1) = S(1);

    %2
    for i = 2:n
        u(i) = A(i,i) - ((A(i,i-1) * A(i-1,i))/(u(i-1)));
        v(i) = S(i) + ((-1*A(i,i-1)*v(i-1))/(u(i-1)));
    end

    %3
    phi = zeros(n,1);
    phi(n) = (v(n))/(u(n));

    %4
    for i = n-1:-1:1
        phi(i) = (1/(u(i))) * (v(i) - (A(i,i+1)*phi(i+1)));
    end
    phi1 = phi;
    
    %Analytical
    C = -13.00868;
    anly_phi = zeros(n,1);
    for i = 1:n
        anly_phi(i) = C*cosh(as(i)/L) + (S0*L^2)/D;
    end
    phi2 = anly_phi;
    
    %Relative Error
    rel_err = zeros(n,1);
    for i = 1:n
        rel_err(i) = abs((phi2(i) - phi1(i))/phi2(i));
    end
    err = [err,max(rel_err)];
end

plot(spaces,err)
set(gca,'FontSize',16)
xlabel('Spaces','FontSize',18)
ylabel('Relative Error','FontSize',18)