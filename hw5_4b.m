n = 5;

A = zeros(n);
b = [];

for i = 1:n
    A(i,i) = 4;
    b(i) = 100;
end
b = b';

for i = 1:(n-1)
    A(i+1,i) = -1;
    A(i,i+1) = -1;
end

start_w = 0;
end_w = 2;
its = 2000;

best_its = Inf;
best_w = [];

iters = [];
ws = linspace(start_w,end_w,its+1);


for i = 1:its+1
    w = ws(i)
    if((w ~= 0) && (w ~= 2))
        [x,its,err] = SOR(A,b,w,[0,0,0,0,0]',10^-6);
        iters(i) = its;
        if(its < best_its)
            best_its = its;
            best_w = [w];
        elseif its == best_its
            best_w(length(best_w)+1) = w;
        end
    else
        iters(i) = NaN;
    end
end

best_w
best_its

plot(ws,iters)
axis([0.95 1.25 0 20])
xlabel('\omega','FontSize',18)
ylabel('Iterations','FontSize',18)
