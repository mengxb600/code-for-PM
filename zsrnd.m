function z=zsrnd(u,a,b,th)
[N,n]=size(u);
mu=ones(N,1)*b+th*a;
sigma=ones(N,n);
y=unifrnd(0,1,N,n);
k=normcdf(0,mu,sigma);
k(k>0.9999)=0.9999;
k(k<0.0001)=0.0001;

x1=k+y.*(1-k);
x2=k.*y;
x=x1.*u+x2.*(1-u);
z=norminv(x,mu,sigma);
end