function I=invtrimatrix(A)
    n=length(A);
    I=diag(1./diag(A));
    for i=1:n
        for j=1:i-1
            disp(I(1:i-1,j))
            disp(A(i,j:i-1).')
            I(i,j)=-dot(I(j:i-1,j),A(i,j:i-1).')/A(i,i);
        end
    end
end
