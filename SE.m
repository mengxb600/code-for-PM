function S1=SE(a,b,c,d,n,nn)
m=eye(n-1);
zo=zeros(n-1,1);
th=linspace(-5.5,5.5,nn);
w=normpdf(th,0,1)*(th(2)-th(1));
for j=1:n
    I=ones(1,n);
    I(j)=0;
    a1=a(I==1);
    b1=b(I==1);
    c1=c(I==1);
    d1=d(I==1);
    p1=c(j)+(d(j)-c(j))*normcdf(a(j)*th+b(j),0,1);
    S=0;
    for k=1:nn
        N=[1,zo'];
        zo1=0*N;
        for i=1:n-1
            p=c1(i)+(d1(i)-c1(i))*normcdf(a1(i)*th(k)+b1(i),0,1);
            q=1-p;
            M=[m*q,zo]+[zo,p*m];
            M=[M;zo1];
            N=N*M;
        end
        S=S+N*(p1(k))*w(k);
    end
    S1(j,:)=S;
end