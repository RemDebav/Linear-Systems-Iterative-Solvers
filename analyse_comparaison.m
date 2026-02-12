ni=20;
n=100;
pas=10;

x=(ni:pas:n);
gi=(ni:pas:n);
j=(ni:pas:n);
sor=(ni:pas:n);
ind=1;
for i=20:pas:n
    A=matrice_definie_positive(i);
    tic
    gauss_implicite(A,rand(i,1),0.001);
    b=toc;
    gi(ind)=b;
    tic
    jacobi(A,rand(i,1),0.001);
    b=toc;
    j(ind)=b;
    tic
    SOR(A,rand(i,1),0.001,0.8);
    b=toc;
    sor(ind)=b;
    ind=ind+1;
end

hold on;
plot(x,gi,'-','DisplayName','Gauss');
plot(x,j,'-','DisplayName','Jacobi')
plot(x,sor,'-','DisplayName','SOR')
hold off;
legend;