function S=SG(a,b,c,d,n,nn)

m=eye(n);
zo=zeros(n,1);
th=linspace(-5.5,5.5,nn);
w=normpdf(th,0,1)*(th(2)-th(1));
S=0;
for k=1:nn
    N=[1,zo'];
    zo1=0*N;
    for i=1:n
        p=c(i)+(d(i)-c(i))*normcdf(a(i)*th(k)+b(i),0,1);
        q=1-p;
        M=[m*q,zo]+[zo,p*m];
        M=[M;zo1];
        N=N*M;
    end
    S=S+N*w(k);
end