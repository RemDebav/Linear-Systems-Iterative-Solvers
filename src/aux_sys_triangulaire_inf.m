function X=systeme_triangulaire_inf(A,b)
    X=zeros(size(b));
    n=length(X);
    for i=1:n
        X(i,1)=(b(i,1)-dot(X(1:i-1,1),A(i,1:i-1)))/A(i,i);
    end
end
