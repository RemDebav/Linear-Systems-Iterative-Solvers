n=50;
sor=linspace(0.1,1.9,n);
x=linspace(0.1,1.9,n);
ind=1;
A=matrice_definie_positive(n);
for omega=linspace(0.1,1.5,n)
    tic
    SOR(A,rand(n,1),0.001,omega);
    b=toc;
    sor(ind)=b;
    x(ind)=omega;
    ind=ind+1;
end
plot(x,sor,'-')

